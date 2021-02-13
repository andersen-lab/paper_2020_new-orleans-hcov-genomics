import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

plt.rcParams['svg.fonttype'] = 'none'

R0 = 2.6
k = 0.16

df = pd.DataFrame(columns = {"run", "generation", "R0", "current_infections"})

ngen = 4

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

df.to_csv("./simulated_cases_4_generations.csv")

f,ax = plt.subplots(figsize=(7.5, 5))
df[df["generation"] == 3]["current_infections"].plot.hist(ax = ax, bins = 500)
ax.axvline(df[df["generation"] == 3]["current_infections"].median(), ls="--", color="red")
ax.axvspan(df[df["generation"] == 3]["current_infections"].quantile(0.025), df[df["generation"] == 3]["current_infections"].quantile(0.975), color="red", alpha = 0.3)
ax.set_title("Total infections after 12 days (4 gen): {} [{} - {}]".format(
    np.round(df[df["generation"] == 3]["current_infections"].median()),
    np.round(df[df["generation"] == 3]["current_infections"].quantile(0.025)),
    np.round(df[df["generation"] == 3]["current_infections"].quantile(0.975))
))
plt.tight_layout()
plt.savefig("simulated_cases.svg", dpi = 300)
plt.close()

# 3rd generation
tmp = df.sort_values("generation").groupby("run").nth(3) - df.sort_values("generation").groupby("run").nth(2)
f,ax = plt.subplots(figsize=(7.5, 5))
ax.hist(tmp["current_infections"], bins = 330)
ax.axvline(tmp["current_infections"].median(), ls="--", color="red")
ax.axvspan(tmp["current_infections"].quantile(0.025), tmp["current_infections"].quantile(0.975), color="red", alpha = 0.3)
ax.set_title("Cases at generation 4: {} [{} - {}]".format(
    np.round(tmp["current_infections"].median()),
    np.round(tmp["current_infections"].quantile(0.025)),
    np.round(tmp["current_infections"].quantile(0.975))
))
plt.tight_layout()
plt.savefig("generation_4.svg", dpi = 300)
plt.close()

f,ax = plt.subplots()
ax.hist(np.random.negative_binomial(k, p, 10000), bins = 100)
plt.tight_layout()
plt.savefig("./secondary_cases.png", dpi=300)
plt.close()


# Nola sum

sum_dist = pd.read_csv("./data/nola_sum_feb_25.csv")

f,ax = plt.subplots()
ax.hist(sum_dist["sum_cases"], bins = 100)
plt.tight_layout()
plt.savefig("./nola_sum_feb_25.png", dpi=300)
plt.close()

# Overlap
infections_simulated_kde = gaussian_kde(df[df["generation"] == 3]["current_infections"].tolist())
infections_sum_kde = gaussian_kde(sum_dist["sum_cases"].tolist())

_x = np.arange(1, sum_dist["sum_cases"].max().astype(int))
_y1 = infections_simulated_kde(_x)
_y2 = infections_sum_kde(_x)
f,ax = plt.subplots(figsize=(7.5, 5))
ax.fill_between(_x, 0, _y1, color="#1874cd", alpha = 0.5)
ax.fill_between(_x, 0, _y2, color="indianred", alpha = 0.5)
ax.axvline(df[df["generation"] == 3]["current_infections"].median(), ls="--", color="#1874cd")
ax.plot(np.minimum(_y1, _y2), color="#000000")
ax.fill_between(_x, 0, np.minimum(_y1, _y2), hatch="///")
ax.axvline(sum_dist["sum_cases"].median(), ls="--", color="indianred")
ax.text(df[df["generation"] == 3]["current_infections"].median(), max(_y1), df[df["generation"] == 3]["current_infections"].median())
ax.text(sum_dist["sum_cases"].median(), max(_y1), sum_dist["sum_cases"].median())
ax.set_title("Overlap: {}".format(np.round(np.trapz(np.minimum(_y1, _y2), _x), 3)))
ax.set_ylim([0, max(_y1)+0.0002])
plt.tight_layout()
plt.savefig("simulated_vs_modelled_dist.svg", dpi = 300)
plt.close()


f,ax = plt.subplots(figsize=(7.5, 5))
df[df["generation"] == 3]["current_infections"].plot.hist(ax = ax, bins = 500, color="#1874cd", alpha = 0.3, density = True)
ax.hist(sum_dist["sum_cases"], bins = 500, color="#6e8b3d", alpha = 0.3, density = True)
ax.axvline(df[df["generation"] == 3]["current_infections"].median(), ls="--", color="red")
ax.axvline(sum_dist["sum_cases"].median(), ls="--", color="#6e8b3d")
ax.set_title("Total infections after 12 days (3 gen): {} [{} - {}]".format(
    np.round(df[df["generation"] == 3]["current_infections"].median()),
    np.round(df[df["generation"] == 3]["current_infections"].quantile(0.025)),
    np.round(df[df["generation"] == 3]["current_infections"].quantile(0.975))
))
plt.tight_layout()
plt.savefig("simulated_vs_modelled_hist.png", dpi =300)
plt.close()

# Doubling time
db = 2
gr = np.log(2)/db

cases = 1 * np.exp(gr * 12)



# 790 - 743 cases

# 13 Feb introduction
# Serial interval of 4 days -> 3 generations

# Get distribution of differences
a = df[df["generation"] == 3].sample(10000, replace=True)
b = sum_dist.sample(10000, replace=True)

diff = b.reset_index(drop=True).values.reshape(1,10000)[0] - a["current_infections"].reset_index(drop = True).values

tmp = diff[(diff >= np.percentile(diff, 2.5)) & (diff <= np.percentile(diff, 97.5))]
kde = gaussian_kde(tmp)

x = np.linspace(min(tmp), max(tmp), 500)
y = kde(x)

f,ax = plt.subplots(figsize=(7.5, 5))
ax.fill_between(x, [0] * len(x), y)
ax.axvline(np.median(diff), ls = "--", color="red")
ax.set_ylim([0, 0.0014])
plt.tight_layout()
plt.savefig("./diff_mardi_gras.svg", dpi = 300)
plt.close()

# COmpute prob
n = [0, 100, 200, 300, 400, 500]
for i in n:
    prob = tmp[tmp > i].shape[0]/tmp.shape[0]
    print("{}: {}%".format(i, np.round(prob * 100, 2)))
