# given 0 <= p1...pn <= 1 with integer sum k, simulate X1...Xn such that in separate each Xi is distributed Bernoulli pi, while combined exactly k are lit

import random
import numpy as np

eps = 0.001

# min i with x <= p_arr[0]+...+p_arr[i]
# not optimized !
def get_interval_index(x, p_arr):
    assert 0 < x and x < sum(p_arr)
    i = -1
    s = 0
    while s < x:
        i += 1
        s += p_arr[i]
    return i

def simulate_cor_bernoulli(p_arr, lit):
    assert 0 <= min(p_arr) and max(p_arr) <= 1 and abs(sum(p_arr) - lit) < eps and lit == int(lit)
    x = random.random()
    return [get_interval_index(x + i, p_arr) for i in range(lit)]

def test(p_arr, lit, rounds = 10**6):
    hist = [0] * len(p_arr)
    for _ in range(rounds):
        sim = simulate_cor_bernoulli(p_arr, lit)
        assert len(sim) == lit and 0 <= min(sim) and max(sim) < len(p_arr)
        for i in range(lit):
            hist[sim[i]] += 1 / rounds
    print(f"p_arr {p_arr}")
    print(f"hist {hist}")
    print(f"test error {max(abs(np.array(hist) - np.array(p_arr)))}")


test([0.5, 0.55, 0.05, 0.2, 0.35, 0.35, 0.33, 0.34, 0.33], 3)
