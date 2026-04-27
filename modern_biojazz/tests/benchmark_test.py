import time
from modern_biojazz.simulation import LocalCatalystEngine
from modern_biojazz.site_graph import ReactionNetwork, Rule, Protein

def create_large_network():
    proteins = {f"P{i}": Protein(name=f"P{i}", sites=[]) for i in range(100)}
    rules = []
    for i in range(1000):
        # some dummy rules
        rules.append(
            Rule(
                name=f"r{i}",
                rule_type="binding",
                reactants=[f"P{i%100}", f"P{(i+1)%100}"],
                products=[f"P{(i+2)%100}"],
                rate=1.0,
            )
        )
    return ReactionNetwork(proteins=proteins, rules=rules)

def test_speed():
    engine = LocalCatalystEngine()
    network = create_large_network()

    start_time = time.time()
    for _ in range(50):
        engine.simulate(network, t_end=1.0, dt=0.01)
    end_time = time.time()
    return end_time - start_time

if __name__ == "__main__":
    t = test_speed()
    print(f"Time: {t:.4f}s")
