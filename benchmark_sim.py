import time
from modern_biojazz.site_graph import ReactionNetwork, Protein, Rule
from modern_biojazz.simulation import LocalCatalystEngine

def create_large_network(n_proteins=2000, n_rules=5000):
    proteins = {f"P{i}": Protein(name=f"P{i}") for i in range(n_proteins)}
    rules = []
    for i in range(n_rules):
        rules.append(Rule(
            name=f"R{i}",
            rule_type="binding",
            reactants=[f"P{i%n_proteins}", f"Token{i}"],
            products=[f"Token{i+1}"],
            rate=1.0
        ))
    return ReactionNetwork(proteins=proteins, rules=rules)

network = create_large_network()

start = time.perf_counter()
# ORIGINAL CODE
species_order = list(network.proteins.keys())
for rule in network.rules:
    for token in [*rule.reactants, *rule.products]:
        if token not in species_order:
            species_order.append(token)

index = {name: i for i, name in enumerate(species_order)}
y0 = [1.0 for _ in species_order]
for i, name in enumerate(species_order):
    if name not in network.proteins:
        y0[i] = 0.0
end = time.perf_counter()
print(f"Old setup time: {end - start:.6f} seconds")

start = time.perf_counter()
# NEW CODE
species_order = list(network.proteins.keys())
seen = set(species_order)
for rule in network.rules:
    for token in [*rule.reactants, *rule.products]:
        if token not in seen:
            seen.add(token)
            species_order.append(token)

index = {name: i for i, name in enumerate(species_order)}
y0 = [1.0 if name in network.proteins else 0.0 for name in species_order]

end = time.perf_counter()
print(f"New setup time: {end - start:.6f} seconds")
