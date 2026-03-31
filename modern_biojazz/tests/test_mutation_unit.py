from __future__ import annotations

import random

from modern_biojazz.mutation import GraphMutator
from modern_biojazz.site_graph import Site


def test_remove_protein_removes_associated_rules(seed_network):
    mutator = GraphMutator(random.Random(1))
    network = seed_network.copy()
    network.rules.append(
        type(network.rules[0])(
            name="bind_temp",
            rule_type="binding",
            reactants=["STAT3", "SOCS3"],
            products=["STAT3:SOCS3"],
            rate=0.1,
        )
    )
    mutator.remove_protein(network, "STAT3")
    assert "STAT3" not in network.proteins
    assert all("STAT3" not in r.reactants for r in network.rules)


def test_action_library_contains_extended_operators(seed_network):
    mutator = GraphMutator(random.Random(2))
    actions = mutator.action_library(seed_network)
    expected = {
        "add_site",
        "add_binding",
        "add_phosphorylation",
        "add_inhibition",
        "remove_rule",
        "modify_rate",
        "remove_site",
        "duplicate_protein",
        "remove_protein",
    }
    assert expected.issubset(set(actions.keys()))


def test_binding_compatibility_requires_both_directions(seed_network):
    mutator = GraphMutator(random.Random(4))
    network = seed_network.copy()
    for site in network.proteins["SOCS3"].sites:
        if site.site_type == "binding":
            site.allowed_partners = ["NOT_STAT3"]
    network.proteins["STAT3"].sites.append(Site(name="b1", site_type="binding", allowed_partners=["SOCS3"]))
    network.proteins["SOCS3"].sites.append(Site(name="b2", site_type="binding", allowed_partners=["NOT_STAT3"]))

    # SOCS3 does not explicitly allow STAT3 here, so compatibility should fail.
    mutator.add_binding_rule(network, "STAT3", "SOCS3")
    assert all(r.rule_type != "binding" for r in network.rules)


def test_remove_protein_cleans_derived_tokens(seed_network):
    mutator = GraphMutator(random.Random(5))
    network = seed_network.copy()
    network.rules.append(
        type(network.rules[0])(
            name="derived_refs",
            rule_type="inhibition",
            reactants=["STAT3_P", "SOCS3"],
            products=["STAT3:SOCS3"],
            rate=0.1,
        )
    )
    mutator.remove_protein(network, "STAT3")
    assert all("STAT3" not in " ".join([*r.reactants, *r.products]) for r in network.rules)

def test_synthesis_and_degradation(seed_network):
    mutator = GraphMutator(random.Random(10))
    network = seed_network.copy()
    mutator.add_synthesis_rule(network, "STAT3")
    mutator.add_degradation_rule(network, "SOCS3")

    syn_rules = [r for r in network.rules if r.rule_type == "synthesis"]
    deg_rules = [r for r in network.rules if r.rule_type == "degradation"]

    assert len(syn_rules) == 1
    assert not syn_rules[0].reactants
    assert syn_rules[0].products == ["STAT3"]

    assert len(deg_rules) == 1
    assert deg_rules[0].reactants == ["SOCS3"]
    assert not deg_rules[0].products

    # Should not raise validation error
    network.validate()

def test_activation_and_dephosphorylation(seed_network):
    mutator = GraphMutator(random.Random(11))
    network = seed_network.copy()
    mutator.add_activation_rule(network, "SOCS3", "STAT3")
    mutator.add_dephosphorylation_rule(network, "SOCS3", "STAT3")

    act_rules = [r for r in network.rules if r.rule_type == "activation"]
    dephos_rules = [r for r in network.rules if r.rule_type == "dephosphorylation"]

    assert len(act_rules) == 1
    assert "STAT3_act" in act_rules[0].products

    assert len(dephos_rules) == 1
    assert "STAT3_P" in dephos_rules[0].reactants
    assert "STAT3" in dephos_rules[0].products

    network.validate()

def test_graph_motifs(seed_network):
    mutator = GraphMutator(random.Random(12))
    network = seed_network.copy()
    mutator.add_protein(network, "P1")
    mutator.add_protein(network, "P2")
    mutator.add_protein(network, "P3")

    mutator.add_feedback_loop(network, "P1", "P2")
    rules_fl = [r for r in network.rules if r.rule_type in ("activation", "inhibition")]
    # Should have P1 activates P2, P2 inhibits P1
    assert any(r.rule_type == "activation" and r.reactants == ["P1", "P2"] for r in rules_fl)
    assert any(r.rule_type == "inhibition" and r.reactants == ["P2", "P1"] for r in rules_fl)

    mutator.add_feedforward_loop(network, "P1", "P2", "P3")
    act_rules = [r for r in network.rules if r.rule_type == "activation"]
    assert any(r.reactants == ["P1", "P3"] for r in act_rules)
    assert any(r.reactants == ["P2", "P3"] for r in act_rules)

    network.validate()
