import gseapy as gp
import pandas as pd
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import numpy as np
import os
import os.path
import mygene
import json


CONTROLS = ["0P", "0M"]
EXPERIMENTALS = ["1A", "2P", "3L", "4C"]

GENESET_DB = "KEGG_2016"

RAW_COUNT_MATRIX_FILE = "combined_counts_new.csv"
COUNT_MATRIX_FILE = "named_combined_counts.csv"
RAW_TEMP = "deseq2_results/{}_raw_ranked_list.csv"
RANKED_TEMP = "deseq2_results/{}_processed_ranked_list.rnk"
OUTDIR = "gsea_results/{}"
BUBBLE_PLOT_TEMP = "bubble_plots/{}_bubble_plot.png"

for direc in ("deseq2_results", "gsea_results", "bubble_plots"):
    if not os.path.isdir(direc):
        os.makedirs(direc)

def rename_count_genes():
    if os.path.exists(COUNT_MATRIX_FILE):
        return
    print("Converting gene IDs for count matrix...")
    raw_matrix_file = RAW_COUNT_MATRIX_FILE
    rnk = pd.read_csv(raw_matrix_file)
    # Convert Ensembl IDs to Gene Symbols
    if True:
        mg = mygene.MyGeneInfo()
        ensembl_ids = rnk["Geneid"].tolist()
        print("  querying mygene...")
        mg_results = mg.querymany(ensembl_ids, scopes="ensembl.gene", fields="symbol", species="human")

        print("  mapping ids...")
        id_map = pd.DataFrame(mg_results)
        # id_map.to_csv("id_map.csv")
    else:
        id_map = pd.read_csv("id_map.csv")
    id_map = id_map[["query", "symbol"]].dropna()
    id_map_dict = id_map.set_index("query")["symbol"].to_dict()
    mapped_rnk = rnk.copy()
    mapped_rnk["Geneid"] = mapped_rnk["Geneid"].map(id_map_dict)
    # mapped_rnk = pd.merge(rnk, id_map, on="Geneid", how="inner")
    mapped_rnk = mapped_rnk.dropna()
    print(f"  Unable to map {len(rnk) - len(mapped_rnk)} of {len(rnk)} gene IDs")
    # convert all gene IDs to upper case
    mapped_rnk["Geneid"] = mapped_rnk["Geneid"].str.upper()
    # remove duplicate gene symbols
    mapped_rnk = mapped_rnk.drop_duplicates(subset="Geneid")

    mapped_rnk.to_csv(COUNT_MATRIX_FILE, index=False)
    print(f"  Successfully converted gene IDs for count matrix")


# if None, no filtering is performed
DESEQ_PVAL_THRESHOLD =[None,]  # [0.05, 0.1, 0.9, None]
def run_dseq2(exp_groups=EXPERIMENTALS, count_threshold=None):
    print("Running DSeq2...")

    # Load the count matrix
    counts = pd.read_csv(COUNT_MATRIX_FILE)
    counts = counts.set_index("Geneid")

    control_cols = [f"{control_group}_col{n}" for control_group in CONTROLS for n in range(1,10)]
    control_cols = [col for col in control_cols if col in counts.columns.tolist()]

    for exp_group in exp_groups:
        print(f"  {exp_group}")
        exp_cols = [f"{exp_group}_col{n}" for n in range(1,10)]
        exp_cols = [col for col in exp_cols if col in counts.columns.tolist()]
        group_cols = control_cols + exp_cols
        
        # filter columns except for exp_col and control columns
        group_counts = counts[group_cols]
        # remove low counts
        olength = group_counts.shape[0]
        if count_threshold is not None:
            group_counts = group_counts[group_counts[group_cols].gt(count_threshold).any(axis=1)]
            print(f"    Removed {olength - len(group_counts)} of {olength} rows with all 0s")

        group_counts = group_counts.T
        print(group_counts)

        metadata = pd.DataFrame({
            "Sample": group_cols,
            "Condition": (["control"] * len(control_cols)) + (["experimental"] * len(exp_cols))
        })
        metadata = metadata.set_index("Sample")
        print(metadata)
        # print(group_counts)

        dds = DeseqDataSet(
            counts=group_counts,
            metadata=metadata,
            design_factors="Condition"
        )
        print(dds)
        dds.deseq2()
        # print(dds)

        stat_res = DeseqStats(dds, contrast=("Condition", "experimental", "control"))
        print(stat_res.summary())

        res = stat_res.results_df
        for threshold in DESEQ_PVAL_THRESHOLD:
            if threshold is None:
                print("  unable to filter")
                break
            osize = res.shape[0]
            res = res[res['pvalue'] < threshold]
            if res.shape[0] < 10000:
                continue
            print(f"  retrieved {osize} results, filtered to {res.shape[0]} with pvalue<{threshold}")

        # save results
        output_file = RAW_TEMP.format(exp_group)
        res.to_csv(output_file)
        print(f"    Results have been saved to {output_file}")


def process_ranked_genes(group):
    raw_list_file = RAW_TEMP.format(group)
    # if not os.path.exists("id_map.csv"):
    print("Loading and processing ranked list...")
    rnk = pd.read_csv(raw_list_file)
    rnk = rnk[["Geneid", "log2FoldChange"]]
    rnk = rnk.dropna()
    # sort by descending log2FoldChange
    rnk = rnk.sort_values(by="log2FoldChange", ascending=False)
    # convert all gene IDs to upper case
    rnk["Geneid"] = rnk["Geneid"].str.upper()

    ranked_list_file = RANKED_TEMP.format(group)
    rnk.to_csv(ranked_list_file, sep="\t", index=False, header=False)
    print(f"  Successfully converted gene IDs for {group}")

def run_gsea(group, ):
    print(f"Running GSEA for {group}...")
    process_ranked_genes(group)

    ranked_list_file = RANKED_TEMP.format(group)
    outdir = OUTDIR.format(group)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    pre_res = gp.prerank(
        rnk=ranked_list_file,
        gene_sets=GENESET_DB,
        threads=4,
        min_size=5,
        max_size=1000,
        permutation_num=1000,
        outdir=outdir,
        seed=6, 
        verbose=True
    )
    print(f"  results saved to {outdir}")


 
# FDR q-value thresholds to try for filtering. If None, no filtering is performed
FDR_THRESHOLDS = [0.05, 0.1, 0.9, None]
# Top N results to include in plot
TOPN = 15
def make_bubble_plot(group, use_log10=True, sortby="NES"):
    def wrap(s, chars=30):
        i = chars
        lines = []
        while i < len(s):
            while i >= 0 and not s[i].isspace():
                i -= 1
            if i == 0:
                i = chars
            lines.append(s[:i].strip())
            s = s[i:].strip()
            i = chars
        return "\n".join(lines)
    
    print(f"Making bubble plot for {group}...")
    gsea_results_file = os.path.join(OUTDIR.format(group), "gseapy.gene_set.prerank.report.csv")
    results = pd.read_csv(gsea_results_file)

    for threshold in FDR_THRESHOLDS:
        if threshold is None:
            print("  unable to filter")
            break
        osize = results.shape[0]
        filtered_results = results[results['FDR q-val'] < threshold]
        if filtered_results.shape[0] < 2:
            continue
        print(f"  retrieved {osize} results, filtered to {results.shape[0]} with FDR<{threshold}")
        results = filtered_results

    # adding stuff for sorting
    if sortby == "NES":
        results["abs_NES"] = results["NES"].abs()
        top = results.sort_values(by="abs_NES", ascending=False).head(TOPN)
    elif sortby == "FDR":
        if use_log10:
            results['-log10_fdr'] = -np.log10(results['FDR q-val'])
            top = results.sort_values(by='', ascending=False).head(TOPN)
        else:
            top = results.sort_values(by='FDR q-val', ascending=True).head(TOPN)
    # top["matched size"] = top['Lead_genes'].apply(lambda x: x.count(";") if isinstance(x, str) else 0)

    plt.figure(figsize=(10, 6))

    if use_log10:
        xs = top['-log10_fdr']
    else:
        xs = top['FDR q-val']
    top["WrappedTerm"] = results["Term"].apply(wrap)
    ys = top["WrappedTerm"]
    dot_sizes = top['Gene %'].str.replace("%", "").astype(float) / 100
    dot_colors = top['NES']  # color by enrichment (Normalized Enrichment Score)

    # min/max normalize dot sizes
    dot_sizes = (dot_sizes - dot_sizes.min()) / (dot_sizes.max() - dot_sizes.min()) + 0.1
    dot_sizes = dot_sizes * 1000

    scatter = plt.scatter(
        x=xs,
        y=ys,
        s=dot_sizes,
        c=dot_colors,
        cmap="coolwarm"
    )
    if use_log10:
        pass
        # plt.gca().invert_xaxis()
    else:
        plt.xlim(0.0, 1.0)

    # NES color bar
    cbar = plt.colorbar(scatter)
    cbar.set_label("Enrichment (NES)", rotation=270, labelpad=15)

    # main GSEA plot
    plt.xlabel("FDR q-value" + (" (-log10)" if use_log10 else ""), fontsize=12)
    plt.title(f"Microglial GSEA, contrast {group} vs controls", fontsize=18)
    plt.grid(axis="both", linestyle="-")
    plt.gca().invert_yaxis()  # invert the y-axis
    plt.tight_layout()
    # plt.margins(x=1.0, y=1.0)

    filename = BUBBLE_PLOT_TEMP.format(group)
    plt.savefig(filename, dpi=300, facecolor="white")
    print(f"  Saved bubble plot to {filename}")

def extract_low_pvals():
    for exp in EXPERIMENTALS:
        fn = f"deseq2_results/{exp}_raw_ranked_list.csv"
        df = pd.read_csv(fn)
        df.set_index("Geneid")
        df = df[(df["padj"] < 0.05) & (df["log2FoldChange"] > 1)]
        df = df.sort_values(by="padj", ascending=True)
        df.to_csv(f"deseq2_results/{exp}_lowest_pvals.csv")

EVAL_FILE = "evaluation_results.json"
def evaluate(params):
    data = {}
    data["params"] = params
    for exp_group in EXPERIMENTALS:
        data[exp_group] = {}
        deseq_results = pd.read_csv(RAW_TEMP.format(exp_group))
        deseq_results["eval_score"] = -np.log10(deseq_results["padj"]) * abs(deseq_results["log2FoldChange"])
        deseq_results = deseq_results.sort_values(by="eval_score", ascending=False)
        # Get the best single eval_score
        
    return

    with open(EVAL_FILE, "w") as f:
        all_data = json.load(f)
        all_data.append(data)
        json.dump(all_data, f, indent=4)
        

if __name__ == "__main__":
    rename_count_genes()
    run_dseq2(count_threshold=1)
    for group in EXPERIMENTALS:
        run_gsea(group)
        make_bubble_plot(group, use_log10=False, sortby="FDR")
