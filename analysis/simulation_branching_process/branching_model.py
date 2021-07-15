import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import itertools
import multiprocessing

plt.rcParams['svg.fonttype'] = 'none'

def simulate_cases(R0, serial_interval):
    print("R0={} SI={}".format(R0, serial_interval))
    k = 0.16
    df = pd.DataFrame(columns = {"run", "generation", "R0", "current_infections"})
    ngen = int(np.round(15/serial_interval))            # Simulate from Feb 11th to Feb25th: 15 days
    for r in range(100000):
        if r % 10000 == 0:
            print(r)
        current_R0 = R0
        population_size = 1270530
        current_infections = 1
        prev_infections = 1
        for i in range(ngen):
            gen_infections = 0
            for l in range(prev_infections):
                p = 1 - (R0 / (R0 + k))
                current_R0 = np.random.negative_binomial(k, p)
                # Ensure each generation has at least 1 transmission and that the first generations always has 1 transmission
                while current_R0 > population_size or (l == prev_infections - 1 and gen_infections < 1 and current_R0 == 0) or (i == 0 and current_R0 == 0):
                    current_R0 = np.random.negative_binomial(k, p)
                    # current_R0 = np.random.poisson(R0)
                    # while current_R0 > population_size:
                    #     current_R0 = np.random.poisson(R0)
                if current_R0 > 0:
                    secondary_infections = current_R0
                    gen_infections += secondary_infections
            prev_infections = gen_infections
            current_infections += gen_infections
            df = df.append({
                "run": r,
                "generation": i,
                "R0": current_R0,
                "current_infections": current_infections
            }, ignore_index=True)
    df.to_csv("./simulated_cases_sr_{}_r0_{}_generations.csv".format(serial_interval, R0))
    df["R0"] = R0
    df["serial_interval"] = serial_interval
    return df

if __name__=="__main__":
    R0s = [2.77, 2.44, 3.17]
    serial_intervals = np.arange(4,9)

    pool = multiprocessing.Pool(processes=15)
    dfs = pool.starmap(simulate_cases, itertools.product(R0s, serial_intervals))
    pool.close()
    pool.join()
    print(len(dfs))

    df = pd.concat(dfs)
    df.to_csv("simulation_100000_results.csv.gz", compression='gzip')



# R0 = 2.6                        # Median R0 = 2.77 (95% HPD: [2.44, 3.17])
# k = 0.16
# serial_interval = 4             # Serial interval from 4 - 8 days.

# df = pd.DataFrame(columns = {"run", "generation", "R0", "current_infections"})

# ngen = np.ceil(15/4)            # Simulate from Feb 11th to Feb25th: 15 days

# for r in range(100000):
#     if r % 10000 == 0:
#         print(r)
#     current_R0 = R0
#     population_size = 1270530
#     current_infections = 1
#     prev_infections = 1
#     for i in range(ngen):
#         gen_infections = 0
#         for l in range(prev_infections):
#             p = 1 - (R0 / (R0 + k))
#             current_R0 = np.random.negative_binomial(k, p)
#             # Ensure each generation has at least 1 transmission and that the first generations always has 1 transmission
#             while current_R0 > population_size or (l == prev_infections - 1 and gen_infections < 1 and current_R0 == 0) or (i == 0 and current_R0 == 0):
#                 current_R0 = np.random.negative_binomial(k, p)
#             # current_R0 = np.random.poisson(R0)
#             # while current_R0 > population_size:
#             #     current_R0 = np.random.poisson(R0)
#             if current_R0 > 0:
#                 secondary_infections = current_R0
#                 gen_infections += secondary_infections
#         prev_infections = gen_infections
#         current_infections += gen_infections
#         df = df.append({
#             "run": r,
#             "generation": i,
#             "R0": current_R0,
#             "current_infections": current_infections
#         }, ignore_index=True)

# df.to_csv("./simulated_cases_sr_{}_r0_{}_generations.csv".format(ngen, serial_interval))
