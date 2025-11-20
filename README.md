# Disparity-in-Differences

The manuscript introducing this method is available at arXiv.

## Getting Started

All core functions and implementations are in the Python script "dispindiffs.py". Especially, this script relies on the Python package "polars" for fast operations on data frames. 

It first needs an edgelist with three columns indicating source, target, and edge weight. 

```python
elist = pl.read_csv(FILEPATH)
Object = dispindiffs.DisparityInDifferences(elist, source=SOURCE_COLUMN_NAME,
      target=TARGET_COLUMN_NAME, weight=WEIGHT_COLUMN_NAME)
Object.calc_disp() # Disparity Filter - Calculate alpha values and add them to a column of elist. 
Object.calc_disp_in_diffs() # Disparity-in-Differences - Calculate alpha values and create a data frame containing calculated values by the method.
```

After these lines were executed, you will have two data frames saved under `Object`. These data frames can be accessed by running `Object.elist` or `Object.elist_did`. 

Now, by using two other methods, you can extract backbones by threshold.
For the disparity filter backbone, use `extr_disp_backbone()`. It returns four outputs: backbone, threshold, number of nodes in the backbone, number of edges in the backbone. For example:
```python
disp_backbone, th, N, E = Object.extr_disp_backbone(th=0.01)    
```
Similarly, for the disparity-in-differences backbone, use `extr_disp_in_diffs_backbone()`. 
```python
disp_in_diffs_backbone, th, N, E = Object.extr_disp_in_diffs_backbone(th=0.01)    
```

All Jupyter notebooks (named "main_eval1_*.ipynb") for evaluation use the same code blocks to extract backbones and then merge external information. 

## Data
Due to the data size limit in Github, files over 25 MB were not uploaded to this repository. Please use the following steps to wrangle data. The Python package "pandas" was mainly used. 

### SciSciNet
Visit https://northwestern-cssi.github.io/sciscinet/ for the original data sets and descriptions. "sciscinet_paperrefs.parquet", "sciscinet_papersources.parquet", and "sciscinet_sources.parquet" were used. A journal citaiton network was constructed through the lines below.
```python
elist = dat.groupby(["citing_paperid", "citing_sourceid", "cited_sourceid"]).size().reset_index(name="c_ij")
elist["f_ij"] = elist["c_ij"] / elist.groupby("citing_paperid")["c_ij"].transform("sum")
elist = elist.groupby(["citing_sourceid", "cited_sourceid"])["f_ij"].sum().reset_index(name="w_ij")
```

### U.S. Airport Network
Download the [2024 DB1B data](https://www.transtats.bts.gov/DL_SelectFields.aspx?gnoyr_VQ=FLM&QO_fu146_anzr=b4vtv0%20n0q%20Qr56v0n6v10%20f748rB). The selected columns are ['ITIN_ID','SEQ_NUM','COUPONS','YEAR','ORIGIN_AIRPORT_ID','QUARTER','ORIGIN','DEST_AIRPORT_ID','DEST','PASSENGERS']. The master coordinate file is available in the same system [LINK](https://www.transtats.bts.gov/DL_SelectFields.aspx?gnoyr_VQ=FLL&QO_fu146_anzr=N8vn6v10%20f722146%20gnoyr5).

### Enron Email Network
Run the included R script "enron_processing_igraphdata.R".

### World Trade Network
The data set "Country Trade by Partner" was downloaded from the [Growth Lab](https://atlas.hks.harvard.edu/data-downloads), Harvard University.
