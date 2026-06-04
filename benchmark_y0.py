import time
import random

species_order = [f"Token{i}" for i in range(1000000)]
network_proteins = {f"Token{i}": True for i in range(500000)}

start = time.perf_counter()
y0 = [1.0 for _ in species_order]
for i, name in enumerate(species_order):
    if name not in network_proteins:
        y0[i] = 0.0
end = time.perf_counter()
print(f"Old approach: {end - start:.6f}s")

start = time.perf_counter()
y0_new = [1.0 if name in network_proteins else 0.0 for name in species_order]
end = time.perf_counter()
print(f"New approach: {end - start:.6f}s")

assert y0 == y0_new
