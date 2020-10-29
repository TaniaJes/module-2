import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# This is the correct formula with the squared std
def welch_t_test(row):
    return (
        (row['ALL_mean'] - row['AML_mean']) / 
        np.sqrt(
            row['ALL_std']**2/len(ALL_SAMP) + row['AML_std']**2/len(AML_SAMP)
        )
    )

df = pd.read_csv(
    'http://pubs.broadinstitute.org/mpr/projects/Leukemia/data_set_ALL_AML_train.txt',
    sep='\t'
    )

cols = df.columns
df = df.reset_index().drop(columns='call.37')
df.columns = cols

df = df.drop(columns=df.columns[df.columns.str.contains('call')])
df = df.drop(columns=['Gene Description'])
df = df.set_index("Gene Accession Number")

fig = plt.figure(figsize=(8, 4), dpi=150)
g = sns.boxenplot(data=df, color='darkgray')
g.set_title("Original Expression")
g.set_xlabel('Sample')
g.set_ylabel('Expression')

fig = plt.figure(figsize=(8, 4), dpi=150)
g = sns.boxenplot(data=df.apply(lambda x: (x - x.mean()) / x.std(), axis=0), color='darkgray')
g.set_title("Sample Standardized Expression")
g.set_xlabel('Sample')
g.set_ylabel('Standardized Expression')

AML_SAMP = [str(x) for x in range(1, 28)]
ALL_SAMP = [str(x) for x in range(28, 39)]

df_std = df.apply(lambda x: (x - x.mean()) / x.std(), axis=0)
df_std[['AML_mean', 'AML_std']] = df_std[AML_SAMP].agg(
    ['mean', 'std'], axis=1
)
df_std[['ALL_mean', 'ALL_std']] = df_std[ALL_SAMP].agg(
    ['mean', 'std'], axis=1
)

df_std['similarity'] = df_std[['ALL_mean', 'ALL_std', 'AML_mean', 'AML_std']].apply(welch_t_test, axis=1)

df_sorted = df_std.sort_values('similarity').drop(columns=['ALL_mean', 'ALL_std', 'AML_mean', 'AML_std', 'similarity'])

fig = plt.figure(figsize=(8, 14), dpi=100)

df_features = df_sorted.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
df_features = df_features.head(25).append(df_features.tail(25))
cmap = sns.diverging_palette(255, 10, s=95, l=35, n=12)

sns.heatmap(
    data=df_features, 
    cmap=cmap, 
    square=1, 
    vmin=-3, 
    vmax=3, 
    cbar_kws={"shrink": .5},
    linewidth=0.1,
    linecolor='gray'
    )