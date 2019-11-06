import iimpute
import pandas as pd


def iimpute_norm(path):
    # read your reads count or RPKM or TPM data
    data = pd.read_csv(path, index_col=0)

    # create I-Impute object
    iimpute_operator = iimpute.IImpute(normalize=False, iteration=True)

    # impute
    imputed_data = iimpute_operator.impute(data)

    return imputed_data
