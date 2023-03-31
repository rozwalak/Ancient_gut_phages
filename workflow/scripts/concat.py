import pandas as pd
import os
from glob import glob


list = glob(snakemake.input.all_csv+'/*.csv')
dfs = [pd.read_csv(f) for f in list]

df_concat = pd.concat(dfs)

df_concat.to_csv(snakemake.output.out_all_csv, header=True, index=False)

cat = ['cat', snakemake.input.all_fasta+'/*.fasta', '>', snakemake.output.out_all_fasta]
os.system(" ".join(cat))
