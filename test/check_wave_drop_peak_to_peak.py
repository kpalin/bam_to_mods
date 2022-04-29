#!/usr/bin/env python

#%%
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

p = Path("/mnt/cg8/home/kpalin/software/bam_to_mods/threshold_tst/wave_data/")


ds = []
for idx, fn in enumerate(
    p.glob(
        "methylation_frequency.Fam_c461_1_19_0711NK_meg.subset.CpG_CTCF_wave.c*.png.tsv"
    )
):
    LIKE_CUTOFF = int(fn.name.split(".")[-3].strip("c"))
    d = pd.read_table(fn)
    d["CUTOFF"] = LIKE_CUTOFF
    ds.append(d)
    del d
    del LIKE_CUTOFF

ds = pd.concat(ds)
# %%
m_mat = ds.set_index(["distance", "CUTOFF"]).percent_meth.unstack()
#%%
mid_desc = m_mat.loc[-20:20].describe()
pre_desc = m_mat.loc[-4000:-2000].describe()
post_desc = m_mat.loc[2000:4000].describe()
#%%

all_desc = pd.DataFrame(
    {"mid": mid_desc.stack(), "pre": pre_desc.stack(), "post": post_desc.stack()}
)

#%%
sns.pairplot(all_desc.loc["mean"])

#%% [markdown]
#
# Pre and post mean methylations are the same, so no need to consider both
#

#%%
sns.jointplot("mid", "pre", hue="CUTOFF", data=all_desc.loc["mean"].reset_index())
# %%
sns.jointplot(
    (all_desc.loc["mean"].reset_index().CUTOFF + 0.5) / 256.0,
    all_desc.loc["mean", "pre"] - all_desc.loc["mean", "mid"],
)
# %%
#%% [markdown]
# 
#
# %%
from IPython.display import Image
Image('../threshold_tst/signal/methylation_frequency.Fam_c461_1_19_0711NK_meg.subset.CpG_CTCF_wave.c125.signal.png')
# %%
Image('../threshold_tst/signal/methylation_frequency.Fam_c461_1_19_0711NK_meg.subset.CpG_CTCF_wave.c253.signal.png')
# %%[markdown]
# From figures above, there seems to be less noise in the 50% (i.e. 125) cutoff image
#%%

D125 = pd.read_table("../threshold_tst/subset_tsv/subset.mC.125.tsv.gz")
# %%
D253 = pd.read_table("../threshold_tst/subset_tsv/subset.mC.253.tsv.gz")

# %%
D125_tot =D125[["uncalled_reads","called_reads"]].sum()
D125_tot*100/D125_tot.sum()
# %%
D253_tot =D253[["uncalled_reads","called_reads"]].sum()
D253_tot*100/D253_tot.sum()

#%% [markdown]
#
# So cutoff 125 will call 83% of reads while 253 will call only 58.6%
#
