from __future__ import annotations

from typing import Dict, Any, List

from modern_biojazz.pipeline import ModernBioJazzPipeline, PipelineConfig
from modern_biojazz.evolution import EvolutionConfig, EvolutionResult
from modern_biojazz.grounding import GroundingResult
from modern_biojazz.site_graph import ReactionNetwork, Rule


class MockEvolutionEngine:
    def __init__(self):
        self.candidate_filter = None
        self.run_called = False
        self.run_kwargs = {}

    def run(self, seed_network: ReactionNetwork, config: EvolutionConfig) -> EvolutionResult:
        self.run_called = True
        self.run_kwargs = {"seed_network": seed_network, "config": config}
        return EvolutionResult(best_network=seed_network, best_score=1.0, history=[0.5, 1.0])


class MockGroundingEngine:
    def __init__(self):
        self.build_constraint_matrix_called = False
        self.match_abstract_to_real_called = False
        self.score_mappings_called = False

    def build_constraint_matrix(self, abstract_types, real_nodes):
        self.build_constraint_matrix_called = True
        return {"STAT3": ["STAT3_HUMAN"]}

    def match_abstract_to_real(self, network, constraints, real_interactions):
        self.match_abstract_to_real_called = True
        return [{"STAT3": "STAT3_HUMAN", "SOCS3": "SOCS3_HUMAN"}]

    def score_mappings(self, mappings, confidence_by_pair):
        self.score_mappings_called = True
        return GroundingResult(mapping=mappings[0], score=0.9, candidates_considered=1)


def test_pipeline_run_without_grounding(seed_network: ReactionNetwork):
    evolution_engine = MockEvolutionEngine()
    grounding_engine = MockGroundingEngine()
    pipeline = ModernBioJazzPipeline(evolution_engine, grounding_engine)

    config = PipelineConfig(
        evolution=EvolutionConfig(generations=2),
        do_grounding=False
    )

    result = pipeline.run(seed_network=seed_network, config=config, grounding_payload=None)

    assert evolution_engine.run_called is True
    assert evolution_engine.candidate_filter is None
    assert result.evolution.best_score == 1.0
    assert result.grounding is None
    assert grounding_engine.build_constraint_matrix_called is False
    assert grounding_engine.match_abstract_to_real_called is False
    assert grounding_engine.score_mappings_called is False


def test_pipeline_run_with_grounding(seed_network: ReactionNetwork, grounding_payload: Dict[str, Any]):
    evolution_engine = MockEvolutionEngine()
    grounding_engine = MockGroundingEngine()
    pipeline = ModernBioJazzPipeline(evolution_engine, grounding_engine)

    config = PipelineConfig(
        evolution=EvolutionConfig(generations=2),
        do_grounding=True
    )

    result = pipeline.run(seed_network=seed_network, config=config, grounding_payload=grounding_payload)

    assert evolution_engine.run_called is True
    assert evolution_engine.candidate_filter is not None
    assert grounding_engine.build_constraint_matrix_called is True
    assert grounding_engine.match_abstract_to_real_called is True
    assert grounding_engine.score_mappings_called is True
    assert result.grounding is not None
    assert result.grounding.score == 0.9


def test_grounding_constraint_filter():
    evolution_engine = MockEvolutionEngine()
    grounding_engine = MockGroundingEngine()
    pipeline = ModernBioJazzPipeline(evolution_engine, grounding_engine)

    allowed = {"A", "B"}
    filter_func = pipeline._grounding_constraint_filter(allowed)

    # Valid network
    network1 = ReactionNetwork()
    network1.rules.append(Rule("r1", "binding", ["A", "B"], ["A:B"], 1.0))
    network1.rules.append(Rule("r2", "phosphorylation", ["A"], ["A_P"], 1.0))
    network1.rules.append(Rule("r3", "inhibition", ["B"], ["B_inh"], 1.0))
    network1.rules.append(Rule("r4", "activation", ["A"], ["A_act"], 1.0))

    assert filter_func(network1) is True

    # Invalid network (unallowed token)
    network2 = ReactionNetwork()
    network2.rules.append(Rule("r1", "binding", ["A", "C"], ["A:C"], 1.0))
    assert filter_func(network2) is False

    # Invalid network (unallowed token suffix)
    network3 = ReactionNetwork()
    network3.rules.append(Rule("r2", "phosphorylation", ["C"], ["C_P"], 1.0))
    assert filter_func(network3) is False
