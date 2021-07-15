from Bio import Phylo
import pandas as pd
from multiprocessing import Pool
import sys
from itertools import repeat

sys.setrecursionlimit(25000)

def get_distance(ids, global_tree, nodes_dict, loc):
    res = [nodes_dict[i] for i in ids if i in nodes_dict]
    anc = global_tree.common_ancestor([i for i in res if i != None])
    dist = 0
    for i in res:
        dist += global_tree.distance(anc, i)
    # # Get paths to MRCA
    # path_mrca = []
    # for i in res:
    #     p = set(anc.get_path(i))
    #     path_mrca.append(p)

    print("{}: {}".format(loc, dist))
    return [dist, len(res)]

if __name__ == '__main__':
    print("Reading JSON ... ")
    full_metadata = pd.read_json("./data/new_api_data.json.gz")

    print("Reading global tree ... ")
    global_tree = Phylo.read("./data/GISAID-hCoV-19-phylogeny-2021-04-02/global.tree", "newick")
    tree_metadata = pd.read_csv("./data/GISAID-hCoV-19-phylogeny-2021-04-02/metadata.csv")

    # Get all nodes and store in dict for quick access
    print("Creating node dict ... ")
    nodes_dict = {}
    for ctr,node in enumerate(global_tree.get_terminals()):
        if ctr % 10000 == 0:
            print(ctr)
        nodes_dict[node.name] = node

    merged_metadata = pd.merge(full_metadata, tree_metadata, left_on="accession_id", right_on="accession_id")

    query = pd.read_csv("./accessin_ids.csv")

    ids = []
    nseqs = []
    grps = []
    grp_col = ["location", "location_id", "division", "division_id", "pangolin_lineage"]
    for name, grp in query.groupby(grp_col):
        ids.append(grp["accession_id"].tolist())
        nseqs.append(grp.shape[0])
        grps.append(name)

    pool = Pool(10)
    tree_res = pool.starmap(get_distance, zip(ids, repeat(global_tree), repeat(nodes_dict), grps))

    pool.close()
    pool.join()

    res = pd.DataFrame(grps, columns = grp_col)
    res["dist"] = [i[0] for i in tree_res]
    res["tree_nseqs"] = [i[1] for i in tree_res]
    res["nseqs"] = nseqs

    res.to_csv("./results_location_lineage.csv")
