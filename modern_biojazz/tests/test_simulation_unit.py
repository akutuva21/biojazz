from __future__ import annotations

import builtins

import urllib.error
from unittest.mock import patch
import pytest

from modern_biojazz.simulation import (
    FitnessEvaluator,
    LocalCatalystEngine,
    UltrasensitiveFitnessEvaluator,
    CatalystHTTPClient,
)


def test_catalyst_http_client_error_path(seed_network):
    client = CatalystHTTPClient(base_url="http://fake", retry_count=2)

    with patch("urllib.request.urlopen") as mock_urlopen, patch("time.sleep") as mock_sleep:
        mock_error = urllib.error.URLError("Connection refused")
        mock_urlopen.side_effect = mock_error

        with pytest.raises(RuntimeError) as exc_info:
            client.simulate(seed_network, t_end=10.0, dt=1.0)

        assert "Failed to simulate network via Catalyst service" in str(exc_info.value)
        # Should attempt 3 times total (1 initial + 2 retries)
        assert mock_urlopen.call_count == 3
        # Should sleep 2 times
        assert mock_sleep.call_count == 2


def test_fitness_evaluator_accepts_backend_network(seed_network):
    engine = LocalCatalystEngine()
    evaluator = FitnessEvaluator(target_output=1.0)
    score = evaluator.score(backend=engine, network=seed_network, t_end=5.0, dt=1.0)
    assert score >= 0.0


def test_ultrasensitive_evaluator_matches_unified_interface(seed_network):
    engine = LocalCatalystEngine()
    evaluator = UltrasensitiveFitnessEvaluator(input_species="STAT3", output_species="SOCS3")
    score = evaluator.score(backend=engine, network=seed_network, t_end=5.0, dt=1.0)
    assert score >= 0.0


def test_local_engine_euler_fallback_when_scipy_unavailable(seed_network, monkeypatch):
    real_import = builtins.__import__

    def fake_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name.startswith("scipy"):
            raise ImportError("forced missing scipy")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    engine = LocalCatalystEngine()
    result = engine.simulate(seed_network, t_end=3.0, dt=1.0, solver="Rodas5P")

    assert result["solver"] == "EulerFallback"
    assert len(result["trajectory"]) == 4
