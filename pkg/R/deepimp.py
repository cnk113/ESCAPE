from deepimpute.deepImpute import MultiNet
import pandas as pd

def deepimp(path)
   data = pd.read_csv(path, index_col=0) # dimension = (cells x genes)
   NN_params = {
        'learning_rate': 1e-4,
        'batch_size': 64,
        'max_epochs': 300,
        'ncores': 5,
        'sub_outputdim': 512,
        'architecture': [
            {"type": "dense", "activation": "relu", "neurons": 200},
            {"type": "dropout", "activation": "dropout", "rate": 0.3}]
    }
    multinet = MultiNet(**NN_params)
    multinet.fit(data,cell_subset=1,minVMR=0.5)
    return multinet.predict(data)
