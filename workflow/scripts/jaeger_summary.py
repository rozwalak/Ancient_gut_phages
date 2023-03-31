import pandas as pd


appended_data = []

for file in glob.glob(snakemake.input):
    data = pd.read_csv(file, sep="\t", header=None)
    appended_data.append(data)

df = pd.concat(appended_data, axis=0, ignore_index=False)
df.to_csv(snakemake.output.table, sep="\t", header=True, index=False)
