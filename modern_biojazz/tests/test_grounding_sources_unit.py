from __future__ import annotations

import json
import urllib.parse
from pathlib import Path
from typing import Any, Dict, List
from unittest.mock import patch, mock_open

from modern_biojazz.grounding_sources import (
    build_grounding_payload_from_sources,
    OmniPathClient,
    INDRAClient,
    load_grounding_snapshot,
)


def test_build_grounding_payload_empty():
    abstract_types: Dict[str, str] = {}
    omnipath_rows: List[Dict[str, Any]] = []
    indra_statements: List[Dict[str, Any]] = []

    payload = build_grounding_payload_from_sources(abstract_types, omnipath_rows, indra_statements)

    assert payload["abstract_types"] == {}
    assert payload["real_nodes"] == []
    assert payload["real_interactions"] == []
    assert payload["confidence_by_pair"] == {}


def test_build_grounding_payload_logic():
    abstract_types = {"STAT3": "tf", "SOCS3": "regulator"}

    # OmniPath rows with different field styles
    omnipath_rows = [
        {"source_genesymbol": "STAT3", "target_genesymbol": "SOCS3_HUMAN"},
        {"source": "STAT3ABC", "target": "SOCS3"},
    ]

    # INDRA statements
    indra_statements = [
        {
            "agents": [{"name": "STAT3_HUMAN"}, {"name": "SOCS3_HUMAN"}],
            "type": "Phosphorylation"
        },
        {
            "agents": [{"name": "STAT3ABC"}, {"name": "SOCS3"}],
            # missing type should default to binding
        }
    ]

    payload = build_grounding_payload_from_sources(abstract_types, omnipath_rows, indra_statements)

    # Verify nodes and type inference
    nodes = {n["name"]: n for n in payload["real_nodes"]}
    assert "STAT3" in nodes
    assert nodes["STAT3"]["type"] == "tf"

    assert "SOCS3_HUMAN" in nodes
    assert nodes["SOCS3_HUMAN"]["type"] == "regulator"

    assert "STAT3ABC" in nodes
    assert nodes["STAT3ABC"]["type"] == "tf"

    assert "SOCS3" in nodes
    assert nodes["SOCS3"]["type"] == "regulator"

    # Verify interactions
    interactions = payload["real_interactions"]
    assert ["STAT3_HUMAN", "SOCS3_HUMAN", "phosphorylation"] in interactions
    assert ["STAT3ABC", "SOCS3", "binding"] in interactions

    # Verify confidence scoring
    conf = payload["confidence_by_pair"]

    # STAT3 (abstract) -> STAT3 (node) : Exact match -> 0.95
    assert conf["STAT3->STAT3"] == 0.95

    # STAT3 (abstract) -> STAT3ABC (node) : Prefix match -> 0.8
    assert conf["STAT3->STAT3ABC"] == 0.8

    # SOCS3 (abstract) -> SOCS3_HUMAN (node) : Prefix with underscore -> 0.95
    assert conf["SOCS3->SOCS3_HUMAN"] == 0.95

    # STAT3 (abstract) -> SOCS3 (node) : No match -> 0.2
    assert conf["STAT3->SOCS3"] == 0.2


@patch("urllib.request.urlopen")
def test_omnipath_client_fetch_interactions(mock_urlopen):
    # Setup mock response
    mock_response = mock_urlopen.return_value.__enter__.return_value
    expected_payload = [{"source": "STAT3", "target": "SOCS3"}]
    mock_response.read.return_value = json.dumps(expected_payload).encode("utf-8")

    client = OmniPathClient()
    genes = ["STAT3", "SOCS3"]
    result = client.fetch_interactions(genes)

    assert result == expected_payload

    # Verify that the URL was constructed properly
    args, kwargs = mock_urlopen.call_args
    req = args[0]
    assert req.method == "GET"

    # URL should contain the genes sorted
    expected_genes = "SOCS3,STAT3"
    assert f"sources={urllib.parse.quote(expected_genes)}" in req.full_url
    assert f"targets={urllib.parse.quote(expected_genes)}" in req.full_url
    assert "genesymbols=1" in req.full_url
    assert req.full_url.startswith("https://omnipathdb.org/interactions/?")


@patch("urllib.request.urlopen")
def test_indra_client_fetch_statements(mock_urlopen):
    mock_response = mock_urlopen.return_value.__enter__.return_value
    expected_statements = [{"type": "Phosphorylation"}]
    mock_response.read.return_value = json.dumps({"statements": expected_statements}).encode("utf-8")

    client = INDRAClient()
    genes = ["STAT3", "SOCS3"]
    result = client.fetch_statements(genes)

    assert result == expected_statements

    args, kwargs = mock_urlopen.call_args
    req = args[0]
    assert req.method == "POST"
    assert req.full_url == "https://api.indra.bio/statements/from_agents"
    assert req.headers["Content-type"] == "application/json"

    payload = json.loads(req.data.decode("utf-8"))
    assert payload["subject"] == genes
    assert payload["object"] == genes
    assert payload["type"] == "Phosphorylation"
    assert payload["format"] == "json"


def test_load_grounding_snapshot():
    expected_data = {"key": "value"}
    json_data = json.dumps(expected_data)

    with patch("builtins.open", mock_open(read_data=json_data)) as mock_file:
        result = load_grounding_snapshot("dummy_path.json")

        assert result == expected_data
        mock_file.assert_called_once_with(Path("dummy_path.json"), "r", encoding="utf-8")
