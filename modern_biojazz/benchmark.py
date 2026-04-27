import time
from modern_biojazz.simulation import LocalCatalystEngine, FitnessEvaluator
from modern_biojazz.site_graph import ReactionNetwork, Protein, Rule

# Create a larger network to amplify the effect
network = ReactionNetwork(metadata={"name": "test"}, proteins={}, rules=[])
for i in range(100):
    network.proteins[f"P{i}"] = Protein(name=f"P{i}", sites=[])
    network.proteins[f"Q{i}"] = Protein(name=f"Q{i}", sites=[])

for i in range(100):
    network.rules.append(Rule(
        name=f"R{i}",
        rule_type="binding",
        rate=1.0,
        reactants=[f"P{i}", f"Q{i}"],
        products=[f"P{i}"]
    ))

engine = LocalCatalystEngine()
evaluator = FitnessEvaluator()

start = time.perf_counter()
# We do 100 runs to get a clearer picture
for _ in range(100):
    engine.simulate(network, t_end=20.0, dt=0.01)
end = time.perf_counter()

print(f"Time: {end - start:.4f}s")
