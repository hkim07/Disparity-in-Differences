# Disparity-in-Differences
A manuscript introducing this method is available at arXiv.

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


