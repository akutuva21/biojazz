from __future__ import annotations

from unittest.mock import MagicMock, patch

from modern_biojazz.benchmarking import benchmark_backend, compare_backends, BenchmarkResult


def test_benchmark_backend():
    mock_backend = MagicMock()
    mock_backend.simulate.return_value = {"trajectory": [{"t": 0, "output": 1.0}]}

    mock_evaluator = MagicMock()
    # We return different scores for different runs to check if the mean is calculated correctly
    mock_evaluator.score.side_effect = [1.0, 2.0, 3.0, 4.0, 5.0]

    mock_network = MagicMock()

    # Mocking time.perf_counter to return specific values
    # Each loop calls perf_counter twice: start, end
    # Run 1: start=0, end=1.0 (duration=1.0)
    # Run 2: start=1.0, end=3.0 (duration=2.0)
    # Run 3: start=3.0, end=6.0 (duration=3.0)
    # Run 4: start=6.0, end=10.0 (duration=4.0)
    # Run 5: start=10.0, end=15.0 (duration=5.0)
    with patch("time.perf_counter", side_effect=[0.0, 1.0, 1.0, 3.0, 3.0, 6.0, 6.0, 10.0, 10.0, 15.0]):
        result = benchmark_backend(
            backend=mock_backend,
            backend_name="TestBackend",
            network=mock_network,
            evaluator=mock_evaluator,
            runs=5,
            t_end=15.0,
            dt=0.5,
            solver="TestSolver",
        )

    assert result.backend_name == "TestBackend"
    assert result.runs == 5
    # mean duration = (1.0 + 2.0 + 3.0 + 4.0 + 5.0) / 5 = 3.0
    assert result.mean_seconds == 3.0
    # mean score = (1.0 + 2.0 + 3.0 + 4.0 + 5.0) / 5 = 3.0
    assert result.mean_score == 3.0

    # Assert simulate was called with correct arguments
    assert mock_backend.simulate.call_count == 5
    mock_backend.simulate.assert_called_with(mock_network, t_end=15.0, dt=0.5, solver="TestSolver")

    # Assert evaluator was called with correct arguments
    assert mock_evaluator.score.call_count == 5
    mock_evaluator.score.assert_called_with(
        simulation_result={"trajectory": [{"t": 0, "output": 1.0}]},
        backend=mock_backend,
        network=mock_network,
        t_end=15.0,
        dt=0.5,
        solver="TestSolver",
    )


def test_benchmark_backend_zero_runs():
    mock_backend = MagicMock()
    mock_evaluator = MagicMock()
    mock_network = MagicMock()

    result = benchmark_backend(
        backend=mock_backend,
        backend_name="ZeroRuns",
        network=mock_network,
        evaluator=mock_evaluator,
        runs=0,
    )

    assert result.backend_name == "ZeroRuns"
    assert result.runs == 0
    assert result.mean_seconds == 0.0
    assert result.mean_score == 0.0
    assert mock_backend.simulate.call_count == 0


@patch("modern_biojazz.benchmarking.benchmark_backend")
def test_compare_backends(mock_benchmark):
    mock_candidate = MagicMock()
    mock_baseline = MagicMock()
    mock_network = MagicMock()
    mock_evaluator = MagicMock()

    # Setup mock returns for benchmark_backend
    mock_benchmark.side_effect = [
        BenchmarkResult(backend_name="candidate", runs=5, mean_seconds=2.0, mean_score=0.8),
        BenchmarkResult(backend_name="baseline", runs=5, mean_seconds=10.0, mean_score=0.9),
    ]

    result = compare_backends(
        candidate=mock_candidate,
        baseline=mock_baseline,
        network=mock_network,
        evaluator=mock_evaluator,
        runs=5,
    )

    # Check that benchmark_backend was called twice with expected args
    assert mock_benchmark.call_count == 2
    mock_benchmark.assert_any_call(mock_candidate, "candidate", mock_network, mock_evaluator, runs=5)
    mock_benchmark.assert_any_call(mock_baseline, "baseline", mock_network, mock_evaluator, runs=5)

    # speedup = baseline_time / candidate_time = 10.0 / 2.0 = 5.0
    assert result == {
        "candidate_mean_seconds": 2.0,
        "baseline_mean_seconds": 10.0,
        "speedup": 5.0,
        "candidate_mean_score": 0.8,
        "baseline_mean_score": 0.9,
    }


@patch("modern_biojazz.benchmarking.benchmark_backend")
def test_compare_backends_zero_candidate_time(mock_benchmark):
    mock_candidate = MagicMock()
    mock_baseline = MagicMock()
    mock_network = MagicMock()
    mock_evaluator = MagicMock()

    # candidate takes 0.0 seconds
    mock_benchmark.side_effect = [
        BenchmarkResult(backend_name="candidate", runs=5, mean_seconds=0.0, mean_score=0.5),
        BenchmarkResult(backend_name="baseline", runs=5, mean_seconds=10.0, mean_score=0.5),
    ]

    result = compare_backends(
        candidate=mock_candidate,
        baseline=mock_baseline,
        network=mock_network,
        evaluator=mock_evaluator,
        runs=5,
    )

    # speedup should be 0.0 to avoid division by zero
    assert result["speedup"] == 0.0
