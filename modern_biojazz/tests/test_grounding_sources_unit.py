from __future__ import annotations

from typing import Any, Dict, List
import json
from unittest.mock import MagicMock, patch
from modern_biojazz.grounding_sources import build_grounding_payload_from_sources, INDRAClient, OmniPathClient


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
def test_indra_client_fetch_statements(mock_urlopen):
    # Setup mock response
    mock_response = MagicMock()
    mock_response.read.return_value = json.dumps({
        "statements": [{"id": 1, "type": "Phosphorylation"}],
        "other_data": "ignored"
    }).encode("utf-8")
    mock_urlopen.return_value.__enter__.return_value = mock_response

    client = INDRAClient()
    genes = ["STAT3", "JAK2"]

    statements = client.fetch_statements(genes, stmt_type="Phosphorylation")

    # Assert correct request structure
    mock_urlopen.assert_called_once()
    req = mock_urlopen.call_args[0][0]

    assert req.full_url == "https://api.indra.bio/statements/from_agents"
    assert req.method == "POST"
    assert req.headers["Content-type"] == "application/json"

    # Assert JSON payload
    payload = json.loads(req.data.decode("utf-8"))
    assert payload == {
        "subject": ["STAT3", "JAK2"],
        "object": ["STAT3", "JAK2"],
        "type": "Phosphorylation",
        "format": "json"
    }

    # Assert returned data matches the statements list
    assert statements == [{"id": 1, "type": "Phosphorylation"}]

@patch("urllib.request.urlopen")
def test_omnipath_client_fetch_interactions(mock_urlopen):
    mock_response = MagicMock()
    mock_payload = [{"source": "A", "target": "B"}]
    mock_response.read.return_value = json.dumps(mock_payload).encode("utf-8")
    mock_urlopen.return_value.__enter__.return_value = mock_response

    client = OmniPathClient()
    genes = ["STAT3", "JAK2", "STAT3"] # Test deduplication and sorting

    interactions = client.fetch_interactions(genes)

    mock_urlopen.assert_called_once()
    req = mock_urlopen.call_args[0][0]

    assert req.method == "GET"
    assert req.full_url.startswith("https://omnipathdb.org/interactions/?")

    # Assert query params
    import urllib.parse
    query = urllib.parse.urlparse(req.full_url).query
    params = urllib.parse.parse_qs(query)

    assert params["genesymbols"] == ["1"]
    assert params["sources"] == ["JAK2,STAT3"] # sorted unique
    assert params["targets"] == ["JAK2,STAT3"]
    assert params["format"] == ["json"]

    assert interactions == mock_payload
