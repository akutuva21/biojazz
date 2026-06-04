from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List
from modern_biojazz.grounding_sources import build_grounding_payload_from_sources, load_grounding_snapshot


def test_load_grounding_snapshot(tmp_path: Path):
    snapshot_data = {"test_key": "test_value", "nodes": [1, 2, 3]}
    file_path = tmp_path / "snapshot.json"

    with open(file_path, "w", encoding="utf-8") as f:
        json.dump(snapshot_data, f)

    # Test with string path
    result_str = load_grounding_snapshot(str(file_path))
    assert result_str == snapshot_data

    # Test with Path object
    result_path = load_grounding_snapshot(file_path)
    assert result_path == snapshot_data


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
