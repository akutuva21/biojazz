import time

def slow_rhs(network, index):
    def rhs(_t: float, y: list[float]) -> list[float]:
        dydt = [0.0 for _ in y]
        for rule in network.rules:
            flux = max(0.0, float(rule.rate))
            for reactant in rule.reactants:
                flux *= max(0.0, y[index[reactant]])
            for reactant in rule.reactants:
                dydt[index[reactant]] -= flux
            for product in rule.products:
                dydt[index[product]] += flux
        return dydt
    return rhs

def fast_rhs_v1(network, index):
    compiled_rules = []
    for rule in network.rules:
        compiled_rules.append((
            max(0.0, float(rule.rate)),
            [index[r] for r in rule.reactants],
            [index[p] for p in rule.products]
        ))

    def rhs(_t: float, y: list[float]) -> list[float]:
        dydt = [0.0] * len(y)
        for rate, r_indices, p_indices in compiled_rules:
            flux = rate
            for r_idx in r_indices:
                flux *= max(0.0, y[r_idx])
            if flux:
                for r_idx in r_indices:
                    dydt[r_idx] -= flux
                for p_idx in p_indices:
                    dydt[p_idx] += flux
        return dydt
    return rhs

def fast_rhs_v2(network, index):
    compiled_rules = []
    for rule in network.rules:
        net_change = {}
        for r in rule.reactants:
            idx = index[r]
            net_change[idx] = net_change.get(idx, 0) - 1
        for p in rule.products:
            idx = index[p]
            net_change[idx] = net_change.get(idx, 0) + 1

        net_change_list = [(idx, float(change)) for idx, change in net_change.items() if change != 0]

        compiled_rules.append((
            max(0.0, float(rule.rate)),
            [index[r] for r in rule.reactants],
            net_change_list
        ))

    def rhs(_t: float, y: list[float]) -> list[float]:
        dydt = [0.0] * len(y)
        for rate, r_indices, net_changes in compiled_rules:
            flux = rate
            for r_idx in r_indices:
                flux *= max(0.0, y[r_idx])
            if flux:
                for idx, change in net_changes:
                    dydt[idx] += flux * change
        return dydt
    return rhs

class MockRule:
    def __init__(self, rate, reactants, products):
        self.rate = rate
        self.reactants = reactants
        self.products = products

class MockNetwork:
    def __init__(self):
        self.rules = []

net = MockNetwork()
index = {}
for i in range(100):
    p, q = f"P{i}", f"Q{i}"
    index[p] = i*2
    index[q] = i*2 + 1
    net.rules.append(MockRule(1.0, [p, q], [p]))

y = [1.0] * 200

rhs1 = slow_rhs(net, index)
s = time.perf_counter()
for _ in range(100000):
    rhs1(0.0, y)
print("Slow:", time.perf_counter() - s)

rhs2 = fast_rhs_v1(net, index)
s = time.perf_counter()
for _ in range(100000):
    rhs2(0.0, y)
print("Fast v1:", time.perf_counter() - s)

rhs3 = fast_rhs_v2(net, index)
s = time.perf_counter()
for _ in range(100000):
    rhs3(0.0, y)
print("Fast v2:", time.perf_counter() - s)
