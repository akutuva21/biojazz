from __future__ import annotations

from unittest.mock import MagicMock, patch

from modern_biojazz.evolution import LLMEvolutionEngine, EvolutionConfig, DeterministicProposer


def test_evaluate_simulation_exception_cegis_feedback(seed_network):
    mock_backend = MagicMock()
    mock_fitness = MagicMock()
    mock_fitness.score.side_effect = Exception("Mock error")

    proposer = DeterministicProposer()

    engine = LLMEvolutionEngine(
        simulation_backend=mock_backend,
        fitness_evaluator=mock_fitness,
        proposer=proposer
    )

    config = EvolutionConfig()

    with patch.object(engine, "_cegis_feedback") as mock_cegis:
        score = engine._evaluate(seed_network, config)

        assert score == 0.0
        mock_cegis.assert_called_once_with(
            network=seed_network,
            score=0.0,
            failure_type="simulation_exception",
            details={"error": "Mock error"}
        )
