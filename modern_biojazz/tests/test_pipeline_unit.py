from __future__ import annotations

import pytest
from unittest.mock import Mock

from modern_biojazz.pipeline import ModernBioJazzPipeline
from modern_biojazz.site_graph import ReactionNetwork, Rule

@pytest.fixture
def pipeline():
    return ModernBioJazzPipeline(evolution_engine=Mock(), grounding_engine=Mock())

def test_grounding_constraint_filter_empty_network(pipeline):
    allowed_symbols = {"A", "B"}
    filter_func = pipeline._grounding_constraint_filter(allowed_symbols)

    network = ReactionNetwork(rules=[])
    assert filter_func(network) is True

def test_grounding_constraint_filter_allowed_exact(pipeline):
    allowed_symbols = {"A", "B"}
    filter_func = pipeline._grounding_constraint_filter(allowed_symbols)

    network = ReactionNetwork(
        rules=[
            Rule(name="r1", rule_type="binding", reactants=["A", "B"], products=["A:B"], rate=1.0)
        ]
    )
    assert filter_func(network) is True

def test_grounding_constraint_filter_allowed_suffixes(pipeline):
    allowed_symbols = {"A", "B"}
    filter_func = pipeline._grounding_constraint_filter(allowed_symbols)

    network = ReactionNetwork(
        rules=[
            Rule(name="r1", rule_type="phosphorylation", reactants=["A"], products=["A_P"], rate=1.0),
            Rule(name="r2", rule_type="inhibition", reactants=["B"], products=["B_inh"], rate=1.0),
            Rule(name="r3", rule_type="activation", reactants=["A"], products=["A_act"], rate=1.0),
        ]
    )
    assert filter_func(network) is True

def test_grounding_constraint_filter_allowed_complex(pipeline):
    allowed_symbols = {"A", "B", "C"}
    filter_func = pipeline._grounding_constraint_filter(allowed_symbols)

    network = ReactionNetwork(
        rules=[
            Rule(name="r1", rule_type="binding", reactants=["A", "B:C"], products=["A:B:C"], rate=1.0)
        ]
    )
    assert filter_func(network) is True

def test_grounding_constraint_filter_disallowed_token(pipeline):
    allowed_symbols = {"A", "B"}
    filter_func = pipeline._grounding_constraint_filter(allowed_symbols)

    network = ReactionNetwork(
        rules=[
            Rule(name="r1", rule_type="binding", reactants=["A", "C"], products=["A:C"], rate=1.0)
        ]
    )
    assert filter_func(network) is False

def test_grounding_constraint_filter_disallowed_suffix(pipeline):
    allowed_symbols = {"A"}
    filter_func = pipeline._grounding_constraint_filter(allowed_symbols)

    network = ReactionNetwork(
        rules=[
            # A_X is not allowed even though A is
            Rule(name="r1", rule_type="modification", reactants=["A"], products=["A_X"], rate=1.0)
        ]
    )
    assert filter_func(network) is False

def test_grounding_constraint_filter_disallowed_complex(pipeline):
    allowed_symbols = {"A", "B"}
    filter_func = pipeline._grounding_constraint_filter(allowed_symbols)

    network = ReactionNetwork(
        rules=[
            Rule(name="r1", rule_type="binding", reactants=["A", "B"], products=["A:C"], rate=1.0)
        ]
    )
    assert filter_func(network) is False
