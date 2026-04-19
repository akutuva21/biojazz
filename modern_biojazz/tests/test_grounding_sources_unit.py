from __future__ import annotations

import json
import urllib.request
from unittest.mock import MagicMock, patch
from typing import Any, Dict, List
from modern_biojazz.grounding_sources import (
    build_grounding_payload_from_sources,
    INDRAClient,
    OmniPathClient,
)


def test_indra_client_fetch_statements():
    client = INDRAClient()
    genes = ["STAT3", "SOCS3"]
    mock_response_data = {"statements": [{"id": "stmt1", "type": "Phosphorylation"}]}

    with patch("urllib.request.urlopen") as mock_urlopen:
        mock_response = MagicMock()
        mock_response.read.return_value = json.dumps(mock_response_data).encode("utf-8")
        mock_response.__enter__.return_value = mock_response
        mock_urlopen.return_value = mock_response

        statements = client.fetch_statements(genes)

        assert statements == mock_response_data["statements"]
        mock_urlopen.assert_called_once()
        req = mock_urlopen.call_args[0][0]
        assert isinstance(req, urllib.request.Request)
        assert req.full_url == "https://api.indra.bio/statements/from_agents"
        assert req.get_method() == "POST"
        assert req.headers["Content-type"] == "application/json"

        sent_payload = json.loads(req.data.decode("utf-8"))
        assert sent_payload["subject"] == genes
        assert sent_payload["object"] == genes
        assert sent_payload["type"] == "Phosphorylation"
        assert sent_payload["format"] == "json"


def test_omnipath_client_fetch_interactions():
    client = OmniPathClient()
    genes = ["STAT3", "SOCS3"]
    mock_response_data = [{"source": "STAT3", "target": "SOCS3"}]

    with patch("urllib.request.urlopen") as mock_urlopen:
        mock_response = MagicMock()
        mock_response.read.return_value = json.dumps(mock_response_data).encode("utf-8")
        mock_response.__enter__.return_value = mock_response
        mock_urlopen.return_value = mock_response

        interactions = client.fetch_interactions(genes)

        assert interactions == mock_response_data
        mock_urlopen.assert_called_once()
        req = mock_urlopen.call_args[0][0]
        assert isinstance(req, urllib.request.Request)
        assert req.get_method() == "GET"
        assert "omnipathdb.org/interactions/" in req.full_url
        assert "genesymbols=1" in req.full_url
        assert "sources=SOCS3%2CSTAT3" in req.full_url  # sorted genes
        assert "targets=SOCS3%2CSTAT3" in req.full_url


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
