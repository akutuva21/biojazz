from __future__ import annotations

import pytest

from modern_biojazz.pipeline import ModernBioJazzPipeline
from modern_biojazz.site_graph import ReactionNetwork, Rule

@pytest.fixture
def pipeline() -> ModernBioJazzPipeline:
    return ModernBioJazzPipeline(None, None)

@pytest.fixture
def dummy_network() -> ReactionNetwork:
    return ReactionNetwork(

        proteins={},
        rules=[],
        metadata={}
    )

def test_filter_exact_match(pipeline, dummy_network):
    dummy_network.rules.append(
        Rule("r1", "binding", reactants=["A", "B"], products=["A:B"], rate=1.0)
    )
    constraint_filter = pipeline._grounding_constraint_filter({"A", "B"})
    assert constraint_filter(dummy_network) is True

def test_filter_exact_match_fails_when_unallowed(pipeline, dummy_network):
    dummy_network.rules.append(
        Rule("r1", "binding", reactants=["A", "C"], products=["A:C"], rate=1.0)
    )
    constraint_filter = pipeline._grounding_constraint_filter({"A", "B"})
    assert constraint_filter(dummy_network) is False

def test_filter_suffixes(pipeline, dummy_network):
    dummy_network.rules.append(
        Rule("r1", "modification", reactants=["A"], products=["A_P"], rate=1.0)
    )
    dummy_network.rules.append(
        Rule("r2", "modification", reactants=["B"], products=["B_inh"], rate=1.0)
    )
    dummy_network.rules.append(
        Rule("r3", "modification", reactants=["C"], products=["C_act"], rate=1.0)
    )
    constraint_filter = pipeline._grounding_constraint_filter({"A", "B", "C"})
    assert constraint_filter(dummy_network) is True

def test_filter_suffixes_fails_when_base_unallowed(pipeline, dummy_network):
    dummy_network.rules.append(
        Rule("r1", "modification", reactants=["D"], products=["D_P"], rate=1.0)
    )
    constraint_filter = pipeline._grounding_constraint_filter({"A", "B", "C"})
    assert constraint_filter(dummy_network) is False

def test_filter_complex_tokens(pipeline, dummy_network):
    dummy_network.rules.append(
        Rule("r1", "binding", reactants=["A", "B"], products=["A:B"], rate=1.0)
    )
    constraint_filter = pipeline._grounding_constraint_filter({"A", "B"})
    assert constraint_filter(dummy_network) is True

def test_filter_complex_tokens_fails_with_suffix(pipeline, dummy_network):
    # According to current implementation, complex tokens split by ":"
    # check if parts are strictly in allowed_symbols. Suffixes like _P inside complex tokens aren't handled properly by the code.
    dummy_network.rules.append(
        Rule("r1", "binding", reactants=["A_P", "B"], products=["A_P:B"], rate=1.0)
    )
    constraint_filter = pipeline._grounding_constraint_filter({"A", "B"})
    assert constraint_filter(dummy_network) is False

def test_filter_complex_tokens_fails_when_part_unallowed(pipeline, dummy_network):
    dummy_network.rules.append(
        Rule("r1", "binding", reactants=["A", "C"], products=["A:C"], rate=1.0)
    )
    constraint_filter = pipeline._grounding_constraint_filter({"A", "B"})
    assert constraint_filter(dummy_network) is False

def test_filter_empty_rules(pipeline, dummy_network):
    constraint_filter = pipeline._grounding_constraint_filter({"A", "B"})
    assert constraint_filter(dummy_network) is True

def test_filter_invalid_suffix(pipeline, dummy_network):
    dummy_network.rules.append(
        Rule("r1", "modification", reactants=["A"], products=["A_unknown"], rate=1.0)
    )
    constraint_filter = pipeline._grounding_constraint_filter({"A", "B"})
    assert constraint_filter(dummy_network) is False

def test_filter_multiple_colons(pipeline, dummy_network):
    dummy_network.rules.append(
        Rule("r1", "binding", reactants=["A", "B:C"], products=["A:B:C"], rate=1.0)
    )
    constraint_filter = pipeline._grounding_constraint_filter({"A", "B", "C"})
    assert constraint_filter(dummy_network) is True
