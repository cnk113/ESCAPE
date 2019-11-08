import doubletdetection
import scrublet as scr
import pandas as pd

def doublets(path, exp=0.1, call=None):
    data = pd.read_csv(path, index_col=0)
    scrub = scr.Scrublet(data, expected_doublet_rate=exp)
    if call:
        predicted_doublets = scr.call_doublets(call)
    else:
        doublet_scores, predicted_doublets = scr.scrub_doublets(min_counts=2, 
                                                                min_cells=3, 
                                                                min_gene_variability_pctl=85, 
                                                                n_prin_comps=30)
    clf = doubletdetection.BoostClassifier(n_iters=50, use_phenograph=False, standard_scaling=True)
    doublets = clf.fit(data).predict(p_thresh=1e-16, voter_thresh=0.5)
    return list(predicted_doublets,doublets)
