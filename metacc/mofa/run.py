from biofam.run.entry_point import entry_point
import pandas as pd
from time import time

infile = "/Users/ricard/data/gastrulation_norsync_stuff/mofa/data_metacc.txt.gz"
outfile = "/Users/ricard/data/gastrulation_norsync_stuff/mofa/hdf5/test.hdf5"

data = pd.read_csv(infile, delimiter="\t", header=0)
lik = ["gaussian"]*data["feature_group"].nunique()

ent = entry_point()
ent.set_data_options(likelihoods=lik, center_features_per_group=True, scale_features=False, scale_views=False)
ent.set_data_df(data)

ent.set_model_options(factors=10, likelihoods=lik, sl_z=False, sl_w=True, ard_z=False, ard_w=True)
ent.set_train_options(iter=1000, convergence_mode="medium", dropR2=None, gpu_mode=True, startELBO=200, elbofreq=5, verbose=False)

ent.build()
ent.run()

ent.save(outfile)