
# Program runs functions from PoolStats.R in python

import os, sys, json
import logging
import pandas as pd
import statistics



def RunPoolStatsPy(mutant_pool_fp, genes_table_fp, nTotalReads, op_dir="."):
    """
    Args:
        mutant_pool_fp (str): Path to mutant pool file
        genes_table_fp (str): Path to genes table
        nTotalReads (int): Number of total reads from all FASTQ files.

    Returns:
        list<success_bool, stats_d>:
            success_bool (bool): True if succeeded, False if failed.
            stats_d (dict): has the following 21 keys: (all 'n' are int, all 'f' are float, 'prcnt' is float)
                'nNonPastEnd', 'nUnique_Insertions', 'nSeenMoreThanOnce', 
                'nSeenExactlyTwice', 'nMedian_strains_per_gene', 'fMean_strains_per_gene', 
                'nMedian_reads_per_gene', 'fMean_reads_per_gene', 
                'fBias_reads_per_gene', 'prcnt_gene_and_transposon_ratio_same_strand', 
                'nCentralInsertions',  'nUniqueCentralInsertions', 
                'nGenesWithInsertions', 'nGenesWithCentralInsertions', 'nGenes',
                'nEssentialGenes', 'nNotEssentialGenes', 'nEssentialGenesHit', 
                'nNotEssentialGenesHit', 'fEssentialGenesHitRatio', 
                'fNonEssentialGenesHitRatio'

    Description:
        We separate pool and genes table by scaffold; for each
        scaffold we check if insertions occured within genes,
        if yes we mark the location and locusId.
        Then we have a new pool dataframe with columns 'f' 
        and 'locusId', which represent insertion fraction
        and name of gene respectively.
        Then we run PoolReport on that dataframe to get
        statistical information.
    """
    
    # We get the two dataframes (pool_df and genes_df) and make sure they
    # have the right columns.
    pool_df = pd.read_table(mutant_pool_fp, sep="\t") 
    if pool_df.shape[0] == 0:
        raise Exception("No rows in mutant pool: " + mutant_pool_fp)
    for col_name in ["nTot", "n", "scaffold", "strand", "pos"]:
        if col_name not in pool_df.columns:
            raise Exception("Mutant pool missing column name: " + col_name)

    genes_df = pd.read_table(genes_table_fp, sep="\t")
    if genes_df.shape[0] == 0:
        raise Exception("No rows in genes table: " + genes_table_fp)
    for col_name in ["scaffoldId", "begin", "end", "strand", "desc"]:
        if col_name not in genes_df.columns:
            raise Exception("Genes Table missing column name: " + col_name)
   
    # Preparing output dict
    final_stats_d = {}
    # We get total number of insertions in the genome
    nonPastEndTrues = pool_df["scaffold"] != "pastEnd"
    nNonPastEnd = nonPastEndTrues.to_list().count(True)
    just_inserts = pool_df[nonPastEndTrues][["scaffold","strand","pos"]]
    # Getting unique inertions
    nUnique_Insertions = just_inserts.drop_duplicates().shape[0]
    info_str = f"{nNonPastEnd} insertions in genome are at {nUnique_Insertions} different locations\n"
    logging.info(info_str)
    nSeenMoreThanOnce = ((nonPastEndTrues) & (pool_df["n"] >= 2)).to_list().count(True)
    nSeenExactlyTwice = ((nonPastEndTrues) & (pool_df["n"] == 2)).to_list().count(True)
    del just_inserts
    # pool_g is like the pool dataframe just with additional columns 'f' and 
    # 'locusId' for insertions that are within genes
    pool_g = create_new_pool_df_with_f_locId(pool_df, genes_df, debug=True) 


    f_list = pool_g["f"].to_list()
    notNA = [x for x in f_list if not pd.isna(x)]
    # If num of rows is 0
    if len(notNA) == 0:
        return [False, {}]

    logging.info("Number of insertions in genes: " + str(len(notNA)))
    pool_report_d = PoolReport(mutant_pool_fp, pool_g, genes_df, nTotalReads, 
                               write_unhit=False, write_surprise=False, write_hitcounts=False, 
                               debug=True)

    final_stats_d = {
        "nNonPastEnd" : nNonPastEnd,
        "nUnique_Insertions" : nUnique_Insertions,
        "nSeenMoreThanOnce" : nSeenMoreThanOnce,
        "nSeenExactlyTwice" : nSeenExactlyTwice
    }
    final_stats_d.update(pool_report_d)

    return [True, final_stats_d]


def PoolReport(mutant_pool_fp, pool_g, genes_df, nTotalReads, 
                write_unhit=False, write_surprise=False, write_hitcounts=False,
               unhit_fp=None, surprise_fp=None, hitcounts_fp=None, 
               debug=False):
    """
    Args:
        pool_df (pandas DataFrame): DataFrame of mutant pool file
        genes_df (pandas DataFrame): DataFrame of genes.GC file
        pool_g (pandas.DataFrame): is like the pool dataframe just with additional 
                               columns 'f' and 
                               'locusId' for insertions that are within genes
        nTotalReads (int): Number of total reads from all FASTQ files.
    Returns:
        pool_report_d (python dict):
            
    """
    # Central 80% of gene
    pool_g = pool_g[pool_g["f"].notna()]
    if debug:
        pool_g.to_csv("debug_pool_g.tsv",sep="\t")

    pool_g2 = pool_g[(pool_g["f"] >= 0.1) & (pool_g["f"] <=0.9)]


    if debug:
        pool_g2.to_csv("debug_pool_g2.tsv",sep="\t")
    # Getting counts of various things:
    nCentralInsertions = pool_g2.shape[0]

    just_central_inserts = pool_g2[["scaffold","strand","pos"]]
    nUniqueCentralInsertions = just_central_inserts.drop_duplicates().shape[0]

    nGenesWithInsertions = len(pool_g["locusId"].unique())

    nGenesWithCentralInsertions = len(pool_g2["locusId"].unique())

    nGenes = genes_df.shape[0]

    genes_df = add_crude_essential_to_genes_df(genes_df)
    genes_df["good_hit"] = [x in pool_g2["locusId"].to_list() for x in genes_df["locusId"].to_list()]
    if debug:
        genes_df.to_csv("debug_updated_genes_df.tsv",sep="\t")

    nEssentialGenes = genes_df["crude_essential"].sum()
    nNotEssentialGenes = (~genes_df["crude_essential"]).sum()
    nEssentialGenesHit = (genes_df["crude_essential"] & genes_df["good_hit"]).sum()
    #essential_locusIds = genes_df[(genes_df["crude_essential"]) & (genes_df["good_hit"])]["locusId"]
    #essential_descs = genes_df[(genes_df["crude_essential"]) & (genes_df["good_hit"])]["desc"]
    nNotEssentialGenesHit = (~genes_df["crude_essential"] & genes_df["good_hit"]).sum()
    fEssentialGenesHitRatio = nEssentialGenesHit/nEssentialGenes
    fNonEssentialGenesHitRatio = nNotEssentialGenesHit/(nGenes - nEssentialGenes)
    logging.debug(f" nEssentialGenesHit: {nEssentialGenesHit},\n"
                  f" nNotEssentialGenesHit: {nNotEssentialGenesHit},\n"
                  f" nEssentialGenes: {nEssentialGenes},\n"
                  f" Total good hit True: {genes_df['good_hit'].sum()},\n"
                  f" nNotEssentialGenes: {nNotEssentialGenes},\n"
                  f" nGenesWithCentralInsertions: {nGenesWithCentralInsertions},\n"
                  f" nGenesWithInsertions: {nGenesWithInsertions},\n"
                  f" fEssential: {fEssentialGenesHitRatio},\n"
                  f" fNonEssential: {fNonEssentialGenesHitRatio}"
                  )
   
    if write_unhit:
        #unhit
        if unhit_fp is None:
            unhit_fp = "unhit_genes.tsv"
        write_unhit_file(genes_df, unhit_fp)

    if write_surprise:
        # surprise
        if surprise_fp is None:
            surprise_fp = "surprise_hits.tsv"
        essential_hits_df = genes_df[(genes_df["crude_essential"]) & (genes_df["good_hit"])]
        essential_hits_df.to_csv(surprise_fp, sep="\t")

    # (HitCounts) Strains & Reads per gene
    locusId_to_strains = pool_g2.groupby("locusId").indices
    locusId_to_strains = {x:len(locusId_to_strains[x]) for x in locusId_to_strains.keys()}
    reads_per_locusId = pool_g2[["locusId", "n"]].groupby("locusId").sum()
    hits_df = pd.DataFrame.from_dict({"locusId": list(reads_per_locusId.index), "nReads": reads_per_locusId["n"],
                "nStrains": [locusId_to_strains[x] for x in reads_per_locusId.index]})
    hits_df.reset_index(drop=True, inplace=True)
    hits_df = hits_df.merge(genes_df, on="locusId")
    hits_df.sort_values(by=["scaffoldId", "begin"], inplace=True)
    if write_hitcounts:
        if hitcounts_fp is None:
            hitcounts_fp = "hitcounts.tsv"
        hits_df.to_csv(hitcounts_fp, sep="\t", index=False)

    median_strains_per_gene = statistics.median(hits_df["nStrains"].to_list())
    mean_strains_per_gene = statistics.mean(hits_df["nStrains"].to_list())
    median_reads_per_gene = statistics.median(hits_df["nReads"].to_list())
    mean_reads_per_gene = statistics.mean(hits_df["nReads"].to_list())
    bias_reads_per_gene = mean_reads_per_gene/median_reads_per_gene
   
    # Transposon and gene on the same strand:
    logging.info("Finding how many transposons and genes are on the same strand.")
    num_same_strand = 0
    genes_df.index = genes_df["locusId"]
    for x in range(pool_g2.shape[0]):
        if genes_df["strand"].loc[pool_g2["locusId"].iloc[x]] == pool_g2["strand"].iloc[x]:
            num_same_strand += 1
    prcnt_gene_and_transposon_ratio_same_strand = (num_same_strand/(pool_g2.shape[0] + 0.001))*100

    
    pool_report_d = {
        "nMedian_strains_per_gene": int(median_strains_per_gene),
        "fMean_strains_per_gene": mean_strains_per_gene,
        "nMedian_reads_per_gene": int(median_reads_per_gene),
        "fMean_reads_per_gene": mean_reads_per_gene,
        "fBias_reads_per_gene": bias_reads_per_gene,
        "prcnt_gene_and_transposon_ratio_same_strand": prcnt_gene_and_transposon_ratio_same_strand,
        "nCentralInsertions": int(nCentralInsertions),
        "nUniqueCentralInsertions": int(nUniqueCentralInsertions),
        "nGenesWithInsertions": int(nGenesWithInsertions),
        "nGenesWithCentralInsertions": int(nGenesWithCentralInsertions),
        "nGenes": int(nGenes),
        "nEssentialGenes": int(nEssentialGenes),
        "nNotEssentialGenes": int(nNotEssentialGenes),
        "nEssentialGenesHit": int(nEssentialGenesHit),
        "nNotEssentialGenesHit": int(nNotEssentialGenesHit),
        "fEssentialGenesHitRatio": fEssentialGenesHitRatio,
        "fNonEssentialGenesHitRatio": fNonEssentialGenesHitRatio
    }


    return pool_report_d

def write_unhit_file(genes_df, op_fp):
    """
    Description:
        Writes out genes with 300 or more nt that had no good hits
        to the file at path op_fp.


	d = without(genes[genes$type==1 & genes$end-genes$begin+1 >= 300 & !genes$goodhit,],words("begin end strand goodhit"));
	d = d[order(d$desc),];
	unhitfile = paste(poolfile,".unhit",sep="");
	writeDelim(d, unhitfile);
	err_printf("Wrote proteins of 300nt or more with no good insertions to %s\n", unhitfile);

    """
    not_good_hits = genes_df[~genes_df["good_hit"]]
    long_enough = not_good_hits[(not_good_hits["end"] - not_good_hits["begin"]) >= 300]
    dropped_labels = long_enough.drop(["begin", "end", "strand", "good_hit"], axis=1)
    dropped_labels.to_csv(op_fp, sep="\t", index=False)
    logging.info("Wrote unhit to " + op_fp)



def add_crude_essential_to_genes_df(genes_df):

    essential_genes_list = []
    for desc in genes_df["desc"].to_list():
        essential = False
        for x in ("ribosomal protein|trna synthetase|cell division").split("|"):
            if x in desc.lower():
                essential = True
                break
        essential_genes_list.append(essential)

    genes_df["crude_essential"] = essential_genes_list
    return genes_df


def create_new_pool_df_with_f_locId(pool_df, genes_df, debug=False):
    """
    Args:
        pool_df (pandas DataFrame): DataFrame of mutant pool file
        genes_df (pandas DataFrame): DataFrame of genes.GC file
    Returns:
        updated_pool_df (None or pd.DataFrame ): The updated pool dataframe
            which now contains 'f' for the fraction of insertion
            if it was inserted into a gene. 'f' will be pd.nan if
            not inserted into a gene.
            It also contains column "locusId" for debugging purposes,
            to check that the insertion was really inserted into
            the gene it is claimed to be inserted into.
    Description:
        We create one new grand pool dataframe composed of 
        the original pool dataframe but split into scaffolds and
        with insertions within genes noted (under columns 'f' and 'locusId').
        For every scaffold we take the subset of the original pool_df,
        take the subset of genes that are in that scaffold as well,
        and then go through them and search for insertions. If an
        insertion is within a gene, then its fraction of insertion
        (what fraction of the gene was the insertion found in, e.g.
        .5 means it was inserted halfway), and note the locusId of
        that gene for debugging purposes.
        We return the grand pool dataframe.
    """

    updated_pool_df = None

    # Both of these are dicts, with scaffold names mapped to 
    # indices of where those scaffolds occured
    pool_scf_group = pool_df.groupby(by="scaffold").indices
    genes_scf_group = genes_df.groupby(by="scaffoldId").indices

    for scf_name in pool_scf_group.keys():
        if scf_name in genes_scf_group:
            logging.info("Running 'add_f_and_locusId_cols' on scaffold " + str(scf_name))
            pool_sub_df = pool_df.iloc[pool_scf_group[scf_name]].copy()
            genes_sub_df = genes_df.iloc[genes_scf_group[scf_name]].copy()
            pool_scf_df = add_f_and_locusId_cols(pool_sub_df, 
                                                 genes_sub_df)
            if debug:
                logging.info("Ran scaffold " + str(scf_name) + " and got nRows:" + str(pool_scf_df.shape[0]))
            if pool_scf_df is not None and not pool_scf_df.shape[0] == 0:
                if updated_pool_df is None:
                    updated_pool_df = pool_scf_df
                else:
                    logging.debug("Appending to updated pool!")
                    # how does this work? Depends on what 'pool_scf_df' is
                    updated_pool_df = updated_pool_df.append(pool_scf_df)

    return updated_pool_df






def add_f_and_locusId_cols(pool_scf_df, genes_scf_df,
                unique=False, minmax=False):
    """
    
    Returns:
        None OR
    
    Description: 
        Takes two dataframes which are the subsets within
        a scaffold of the pool dataframe and the genes dataframe.
        Returns a dataframe which
        looks at insertions within a given scaffold in
        pool_scf_df and finds if it occured within genes_scf_df.
        'pos' from pool_scf_df must be between 'begin'/'end' in
        genes_scf_df
    TD:
        Replace 'x_df' with pool_scf_df and
        'y_df' with genes_scf_df
    """
    if pool_scf_df.shape[0] == 0 or genes_scf_df.shape[0] == 0:
        return None

    pool_ix = 0
    gene_ix = 0

    pool_scf_df.sort_values(by="pos", ascending=True, inplace=True)
    genes_scf_df.sort_values(by="begin", ascending=True, inplace=True)
    nPool_scf = pool_scf_df.shape[0]
    
    # This list stores the insert frac for each insert (pd.NA if not in a gene) 
    insert_frac_list = []
    # This list stores the locusId of the inserted gene for debugging
    locusIds = []
    while gene_ix < genes_scf_df.shape[0] and pool_ix < nPool_scf:
        crt_gene = genes_scf_df.iloc[gene_ix] 
        crt_insert = pool_scf_df.iloc[pool_ix]
        if crt_insert["pos"] <= crt_gene["begin"]:
            pool_ix += 1
            insert_frac_list.append(pd.NA)
            locusIds.append(pd.NA)
        elif crt_insert["pos"] >= crt_gene["end"]:
            gene_ix += 1
        else:
            fInsert = (crt_insert["pos"] - crt_gene["begin"]) / (crt_gene["end"] - crt_gene["begin"])
            insert_frac_list.append(fInsert)
            locusIds.append(crt_gene["locusId"])
            pool_ix += 1

    
    if pool_ix < nPool_scf:
        insert_frac_list += [pd.NA]*(nPool_scf - pool_ix)
        locusIds += [pd.NA]*(nPool_scf - pool_ix)

    pool_scf_df["f"] = insert_frac_list
    pool_scf_df["locusId"] = locusIds 

    return pool_scf_df
            

     
def PSpause(val):
    raise Exception(str(val))



def main():
    args = sys.argv
    if args[-1] != "1":
        help_str = "python3 PoolStats.py pool genes.tab totalReads op_dir 1\n" + \
	            "   If successful, writes metrics to op_dir/pool.stats and" + \
                    "   writes the set of unhit protein-coding genes (no good " + \
                    "insertion) to op_dir/pool.unhit."
        print(help_str)
        sys.exit(1)
    else:
        pool_fp = args[1]
        genes_fp = args[2]
        nTotalReads = int(args[3])
        op_dir = args[4]
        logging.basicConfig(level=logging.DEBUG)
        res = RunPoolStatsPy(pool_fp, genes_fp, nTotalReads, op_dir=op_dir)
        print(res[0])
        print(json.dumps(res[1], indent=2))
        print(res[1].keys())
        
if __name__ == "__main__":
    main()
