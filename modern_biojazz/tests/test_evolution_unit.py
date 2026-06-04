from __future__ import annotations

import random
from typing import Any, Dict

import pytest

from modern_biojazz.evolution import (
    DeterministicProposer,
    EvolutionConfig,
    EvolutionResult,
    LLMEvolutionEngine,
)
from modern_biojazz.mutation import GraphMutator
from modern_biojazz.site_graph import ReactionNetwork


class MockBackend:
    def simulate(
        self,
        network: ReactionNetwork,
        t_end: float,
        dt: float,
        solver: str,
        initial_conditions: Dict[str, float] | None = None,
    ) -> Dict[str, Any]:
        return {"trajectory": [{"t": 0.0, "output": 0.0}, {"t": t_end, "output": 1.0}]}


class MockFitness:
    def __init__(self):
        self.evals = 0

    def score(
        self,
        simulation_result: Dict[str, Any] | None = None,
        *,
        backend: Any = None,
        network: ReactionNetwork | None = None,
        t_end: float = 20.0,
        dt: float = 1.0,
        solver: str = "Rodas5P",
        initial_conditions: Dict[str, float] | None = None,
    ) -> float:
        self.evals += 1
        return 0.1 + (self.evals * 0.001)


def test_evolution_engine_basic_run(seed_network):
    backend = MockBackend()
    fitness = MockFitness()
    proposer = DeterministicProposer(feedback_log=[])
    mutator = GraphMutator()
    rng = random.Random(42)

    engine = LLMEvolutionEngine(
        simulation_backend=backend,
        fitness_evaluator=fitness,
        proposer=proposer,
        mutator=mutator,
        rng=rng,
    )

    config = EvolutionConfig(
        population_size=4,
        generations=3,
        mutations_per_candidate=2,
        islands=2,
        migration_interval=2,
        migration_count=1,
        sim_t_end=1.0,
        sim_dt=1.0,
    )

    result = engine.run(seed_network, config)

    assert isinstance(result, EvolutionResult)
    assert len(result.history) == 4  # Initial + 3 generations
    assert result.history[0] <= result.history[-1]
    assert fitness.evals > 0
    # The proposer feedback should be recorded at the end of each generation
    assert len(proposer.feedback_log) >= 3


class FailingFitness:
    def score(self, *args, **kwargs) -> float:
        raise RuntimeError("Simulation crashed")


def test_evolution_engine_simulation_exception(seed_network):
    fitness = FailingFitness()
    proposer = DeterministicProposer(feedback_log=[])
    engine = LLMEvolutionEngine(
        simulation_backend=MockBackend(),
        fitness_evaluator=fitness,
        proposer=proposer,
    )

    config = EvolutionConfig(population_size=2, generations=1, islands=1)
    result = engine.run(seed_network, config)

    assert result.best_score == 0.0
    assert any("simulation_exception" in log for log in proposer.feedback_log)


def test_evolution_engine_candidate_filter(seed_network):
    def reject_all(network: ReactionNetwork) -> bool:
        return False

    proposer = DeterministicProposer(feedback_log=[])
    engine = LLMEvolutionEngine(
        simulation_backend=MockBackend(),
        fitness_evaluator=MockFitness(),
        proposer=proposer,
        candidate_filter=reject_all,
    )

    config = EvolutionConfig(population_size=2, generations=1, islands=1)
    engine.run(seed_network, config)

    assert engine.filter_rejection_count > 0
    assert any("filter_rejection" in log for log in proposer.feedback_log)


class ZeroFitness:
    def score(self, *args, **kwargs) -> float:
        return 0.0


def test_evolution_engine_zero_fitness_feedback(seed_network):
    fitness = ZeroFitness()
    proposer = DeterministicProposer(feedback_log=[])
    engine = LLMEvolutionEngine(
        simulation_backend=MockBackend(),
        fitness_evaluator=fitness,
        proposer=proposer,
    )

    config = EvolutionConfig(population_size=2, generations=1, islands=1)
    engine.run(seed_network, config)

    assert any("low_fitness" in log for log in proposer.feedback_log)


class NoFeedbackProposer:
    def propose(self, model_code: str, action_names: list[str], budget: int) -> list[str]:
        return []


def test_evolution_engine_no_feedback_proposer(seed_network):
    proposer = NoFeedbackProposer()
    engine = LLMEvolutionEngine(
        simulation_backend=MockBackend(),
        fitness_evaluator=MockFitness(),
        proposer=proposer,
    )

    config = EvolutionConfig(population_size=2, generations=1, islands=1)
    # This should not raise an exception
    engine.run(seed_network, config)
