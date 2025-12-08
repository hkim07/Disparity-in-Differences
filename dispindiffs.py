import polars as pl
import numpy as np
from scipy.stats import beta
from tqdm import tqdm

class DisparityInDifferences:
    def __init__(self, elist, source="source", target="target", weight="weight", n_samples=10000):
        """
        elist: an edgelist in the Polars dataframe format.
        source: the column name for sources
        target: the column name for targets
        weight: the column name for edge weights
        """
        self.elist = elist
        self.source = source
        self.target = target
        self.weight = weight
        self.n_samples = n_samples
        self.k_out = None
        self.elist_did = None

    def preprocessing(self):
        # Rename source, target, and weight columns to "source", "target", and "weight"
        self.elist = self.elist.rename({self.source: "source", self.target: "target", self.weight: "weight"})
        
        # Remove all self-loops and obtain a dictionary for out-degrees        
        self.elist = self.elist.with_columns(
            pl.col("source").cast(pl.Utf8),
            pl.col("target").cast(pl.Utf8)
        )
        self.elist = self.elist.filter(pl.col("source") != pl.col("target"))

        self.k_out = self.elist.group_by("source").agg(pl.col("target").n_unique().alias("k_out"))
        self.k_out = dict(zip(self.k_out["source"], self.k_out["k_out"]))
        
        # Normalize edge weights to calculate p_ij
        sum_w_ij = (
            self.elist.group_by("source")
                 .agg(pl.col("weight").sum().alias("sum_w_ij"))
        )
        self.elist = self.elist.join(sum_w_ij, on="source")
        self.elist = self.elist.with_columns((pl.col("weight") / pl.col("sum_w_ij")).alias("p_ij"))
        
    def calc_disp(self):
        if "p_ij" not in self.elist.columns:
            self.preprocessing()
        self.elist = self.elist.with_columns(
            pl.col("source").replace_strict(self.k_out, default=None).alias("k_i_out")
        )
        self.elist = self.elist.with_columns(
            ((1 - pl.col("p_ij")) ** (pl.col("k_i_out") - 1)).alias("disp_alpha"),
        )
        
    def extr_disp_backbone(self, th=0.05):
        disparity_backbone = self.elist.filter(pl.col("disp_alpha")<=th)
        N = len(set(disparity_backbone["source"].unique()) | set(disparity_backbone["target"].unique()))
        E = len(disparity_backbone)
        return disparity_backbone, th, N, E

    def calc_disp_in_diffs(self):
        if "p_ij" not in self.elist.columns:
            self.preprocessing()
        print("Merging bilateral relations")
        self.elist_did = self.elist.with_columns(
            (pl.col("source")+"/"+pl.col("target")).alias("original_key"),
            pl.when(pl.col("source") < pl.col("target"))
                .then(pl.col("source")+"/"+pl.col("target"))
                .otherwise(pl.col("target")+"/"+pl.col("source"))
                .alias("sorted_key")
        )

        self.elist_did = self.elist_did.with_columns(
            (pl.col("original_key")==pl.col("sorted_key")).alias("key_matched")
        )
        self.elist_did = self.elist_did.sort(["sorted_key", "key_matched"], descending=True)
        
        self.elist_did = self.elist_did.group_by("sorted_key").agg([
            pl.len().alias("n"),
            pl.col("source").first().alias("i"),
            pl.col("target").first().alias("j"),
            pl.col("p_ij").first().alias("p_ij"),
            pl.col("p_ij").last().alias("p_ji"),
        ])
        
        self.elist_did = self.elist_did.filter(pl.col("n")==2)        
        self.elist_did = self.elist_did.with_columns([
            (pl.col("p_ij") - pl.col("p_ji")).alias("D_ij"),
            pl.col("i").replace_strict(self.k_out, default=0).alias("k_i_out"),
            pl.col("j").replace_strict(self.k_out, default=0).alias("k_j_out")
        ])        

        self.elist_did = self.elist_did.sort("sorted_key").drop("sorted_key")    
        print("Generating pre-sampled values from beta distributions")
        uniq_k_out_values = set(self.elist_did["k_i_out"].unique()) | set(self.elist_did["k_j_out"].unique())
        pre_samples = {}
        for k in uniq_k_out_values:
            try:
                pre_samples[k] = beta.rvs(1, k-1, size=self.n_samples, random_state=k)
            except:
                continue
                
        print("Calculating statistical significance")
        sig = []
        for row in tqdm(self.elist_did.iter_rows(named=True), total=len(self.elist_did)):
            try:
                E_ij = pre_samples[row["k_i_out"]] - pre_samples[row["k_j_out"]]
                alpha = np.mean(row["D_ij"] > E_ij)
                sig.append(alpha)
            except:
                sig.append(None)
        self.elist_did = self.elist_did.with_columns(pl.Series("disp_in_diffs_alpha", sig))        
        print("Done")

    def extr_disp_in_diffs_backbone(self, th=0.05):
        """
        th/2 will be applied to each side
        """
        tmp = self.elist_did.filter((pl.col("disp_in_diffs_alpha")<=(th/2)) | (pl.col("disp_in_diffs_alpha")>=(1-th/2)))
        did_backbone = []
        for row in tmp.iter_rows(named=True):
            if row["D_ij"]<0:
                did_backbone.append((row["j"], row["i"]))
            if row["D_ij"]>=0:
                did_backbone.append((row["i"], row["j"]))     
        did_backbone = pl.DataFrame(did_backbone, schema=["source", "target"], orient="row")    
        N = len(set(did_backbone["source"].unique()) | set(did_backbone["target"].unique()))
        E = len(did_backbone)
        return did_backbone, th, N, E