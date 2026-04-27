import pytest
from unittest.mock import MagicMock

from modern_biojazz.pipeline import ModernBioJazzPipeline
from modern_biojazz.site_graph import ReactionNetwork, Rule


@pytest.fixture
def pipeline():
    evolution_engine_mock = MagicMock()
    grounding_engine_mock = MagicMock()
    return ModernBioJazzPipeline(evolution_engine_mock, grounding_engine_mock)


def test_grounding_constraint_filter_is_allowed_token(pipeline):
    allowed_symbols = {"A", "B", "C"}
    constraint_filter = pipeline._grounding_constraint_filter(allowed_symbols)

    # We need to access the inner _is_allowed_token function
    # The returned value is _filter(network), which captures _is_allowed_token in its closure.
    # To test the tokens directly, we can test _filter with dummy networks containing one token.

    def check_token(token):
        network = ReactionNetwork(
            proteins={},
            rules=[Rule(name="r1", rule_type="binding", reactants=[token], products=[], rate=1.0)]
        )
        return constraint_filter(network)

    # Allowed tokens
    assert check_token("A") is True
    assert check_token("B") is True

    # Suffixed tokens
    assert check_token("A_P") is True
    assert check_token("B_inh") is True
    assert check_token("C_act") is True

    # Complex tokens
    assert check_token("A:B") is True
    assert check_token("A:B:C") is True
    assert check_token("A_P:B") is False # The function logic splits by : and checks exact parts, so A_P won't be allowed as part of A:B unless A_P is in allowed
    assert check_token("A:B_P") is False # same reason

    # Disallowed tokens
    assert check_token("D") is False
    assert check_token("A_X") is False
    assert check_token("D_P") is False
    assert check_token("D_inh") is False
    assert check_token("D_act") is False
    assert check_token("A:D") is False


def test_grounding_constraint_filter_filter(pipeline):
    allowed_symbols = {"A", "B"}
    constraint_filter = pipeline._grounding_constraint_filter(allowed_symbols)

    network_allowed = ReactionNetwork(
        proteins={},
        rules=[
            Rule(name="r1", rule_type="binding", reactants=["A", "B"], products=["A:B"], rate=1.0),
            Rule(name="r2", rule_type="phosphorylation", reactants=["A"], products=["A_P"], rate=1.0),
            Rule(name="r3", rule_type="inhibition", reactants=["B"], products=["B_inh"], rate=1.0)
        ]
    )
    assert constraint_filter(network_allowed) is True

    network_disallowed = ReactionNetwork(
        proteins={},
        rules=[
            Rule(name="r1", rule_type="binding", reactants=["A", "C"], products=["A:C"], rate=1.0),
        ]
    )
    assert constraint_filter(network_disallowed) is False

    network_partially_disallowed = ReactionNetwork(
        proteins={},
        rules=[
            Rule(name="r1", rule_type="binding", reactants=["A", "B"], products=["A:B"], rate=1.0),
            Rule(name="r2", rule_type="phosphorylation", reactants=["C"], products=["C_P"], rate=1.0),
        ]
    )
    assert constraint_filter(network_partially_disallowed) is False
