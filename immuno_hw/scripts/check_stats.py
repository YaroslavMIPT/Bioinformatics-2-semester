import os
import pandas as pd
from repseq import clonosets as cl
from repseq import stats

ucb_dir = "data/ucb_figshare"
clonosets = cl.find_all_mixcr_clonosets(ucb_dir, show_offtarget=False, offtarget_threshold=0.01)
metadata = pd.read_csv(os.path.join(ucb_dir, "ucb_figshare_TRB_metadata.csv"))[
    ["sample_id", "experimental_group", "donor_id", "subset"]
]

for sub in ["nCD4", "nCD8"]:
    print(f"\n================ {sub} ================")
    df_sub = clonosets.merge(metadata).query("subset == @sub").reset_index(drop=True)
    print("Образцов по группам:")
    print(df_sub["experimental_group"].value_counts())
    st = stats.calc_clonoset_stats(df_sub, cpu=4).merge(metadata[["sample_id", "experimental_group"]])
    print("\n10 самых маленьких образцов по reads_func:")
    print(st[["sample_id", "experimental_group", "clones_func", "reads_func"]].sort_values("reads_func").head(10).to_string(index=False))

print("\n================ BCR (MiXCR) ================")
bcr_cl = cl.find_all_mixcr_clonosets("outputs/bcr_mixcr", show_offtarget=False, offtarget_threshold=0.01)
print(bcr_cl)
if len(bcr_cl) > 0:
    bcr_st = stats.calc_clonoset_stats(bcr_cl, cpu=2)
    print(bcr_st[["sample_id", "chain", "clones", "clones_func", "reads", "reads_func"]].to_string(index=False))
