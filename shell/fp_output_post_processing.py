import pandas as pd
import seaborn as sns
from pathlib import Path
import numpy as np
import argparse
import time
import matplotlib.pyplot as plt

from contaminants import CONTAMINANTS
#Utility functions and classes for MS data import

def setClass(accession, contaminants_list=CONTAMINANTS):
    '''
    assign a class to the peptide
    NEOTRANSCRIPT is the tag that is added to neotranscript entries by the add_PE.py script
    lcl is a tag added to the result of translation with ORFfinder
    sp is the tag of proteins from swissprot (uniprot)
    source of the contaminants list: https://www.thegpm.org/crap/
    '''
    if ((('NEOTRANSCRIPT' in accession) or ('lcl' in accession)) and 'sp|' not in accession):
        return("Non canonical")
    elif ((('NEOTRANSCRIPT' in accession) or ('lcl' in accession)) and 'sp|' in accession):
        return("Ambiguous")
    elif any(contaminant in accession for contaminant in contaminants_list):
        return('Contaminant')
    else:
        return('Canonical')


class MSrun:
    def __init__(self, metadata_row):
        self.sample = metadata_row["sample"]
        self.DB = metadata_row["DB"]
        self.search_engine = metadata_row["search_engine"]
        self.fdr = metadata_row["fdr"]
        self.hla = metadata_row["hla_alleles"].split("|")

        search_engine = metadata_row["search_engine"]
        #print(metadata_row["path"])
        if "Fragpipe" in search_engine:
            self.psm = self.from_fragpipe(metadata_row["path"])
        else:
            raise ValueError(f"Provided search engine ({search_engine}) is not supported.\n \
            Only Fragpipe search engine is supported.")
        self.unique = self.psm.drop_duplicates(subset=['seq_clear'], inplace=False)

    @staticmethod
    def from_fragpipe(path):
        def get_fraction(row):
            return row["Spectrum"].split(".")[0]
	
        def get_scan(row):
            return int(row["Spectrum"].split(".")[2])

        def extract_accessions(row):
            if not (row["Is Unique"] | isinstance(row["Mapped Proteins"], float)):
                accessions = ",".join((row["Protein"],row["Mapped Proteins"]))
                accessions = accessions.replace(",",";")
                accessions = accessions.replace(" ","")
            else:
                accessions = row["Protein"]
            return str(accessions)

        psm_file = pd.read_csv(path, sep="\t")
        #print(len(psm_file.index))
        dict_info = {"type": ["PEPTIDE"]*len(psm_file),
            "rt": psm_file["Retention"], # convert into seconds
            #"expectation": psm_file["Expectation"], # see http://pappso.inrae.fr/bioinfo/i2masschroq/documentation/html/chap_fundamentals-bottom-up-proteomics.html E-value
            #"spectral_entropy": psm_file["SpectralSim"], # another quality metric for spectra see https://doi.org/10.1038/s41592-021-01331-z 
            "score": psm_file["Hyperscore"], # using PEP score as q-value
            "rank": [None]*len(psm_file),
            "sequence": psm_file["Modified Peptide"],
            "seq_clear": psm_file["Peptide"],
            "aa_before": psm_file["Prev AA"],
            "aa_after": psm_file["Next AA"],
            "score_type": ["PSM"]*len(psm_file),
            "search_identifier": ["Fragpipe"]*len(psm_file),
            "accessions": psm_file.apply(extract_accessions, axis=1),
            "start": psm_file["Protein Start"],
            "end": psm_file["Protein End"],
            "fraction": psm_file.apply(get_fraction, axis=1),
            "scan": psm_file.apply(get_scan, axis=1),
            #"intensities": psm_file["Intensities"]
            #    "sample": meta["sample"],
            #"fdr": [None]*len(psm_file)
            }
        psm = pd.DataFrame.from_dict(dict_info)
        return psm

    def isCanonical(self):
        self.psm["class"] = self.psm.accessions.apply(lambda x: setClass(x))
        self.unique["class"] = self.unique.accessions.apply(lambda x: setClass(x))

#Utility functions for plotting

def combine(counts, percentages):
    fmt = "{}\n{:.1f}%".format
    return [fmt(c, p) for c, p in zip(counts, percentages)]

def plot_class_distribution(run_, max_, save_dir=None):
    data_ = run_.unique.copy() 
    data_['class'] = pd.Categorical(data_["class"], categories=["Canonical", "Ambiguous", "Non canonical", "Contaminant"])

    plt.figure(figsize=(6, 3))
    ax = sns.countplot(data=data_, x='class', palette={"Canonical": "#bcbcbc", "Ambiguous": "#b63fc0", "Non canonical": "#eaaa6a", "Contaminant": "#9e6156"}, width=0.75)
    counts = data_['class'].value_counts()[["Canonical", "Ambiguous", "Non canonical", "Contaminant"]]
    percentages = (counts / counts.sum()) * 100
    annotations = combine(counts.values, percentages.values)

    # ax.patches only has bars for present categories — match by x-position to category index
    categories = ["Canonical", "Ambiguous", "Non canonical", "Contaminant"]
    for p in ax.patches:
        cat_index = round(p.get_x() + p.get_width() / 2)  # x tick position = category index
        annotation = annotations[cat_index]
        ax.annotate(annotation, (p.get_x() + p.get_width() / 2., p.get_height()),
                    ha='center', va='bottom', fontsize=9)

    ax.set_ylim(0, max_)
    ax.set_ylabel("Count")
    ax.set_title(f"{run_.sample} {run_.search_engine} peptides count FDR = {run_.fdr}")
    plt.xticks(rotation=0)

    ax.yaxis.grid(True, linestyle='--', linewidth=0.5, color='lightgray')
    
    sns.despine()

    if save_dir:
        plt.savefig(Path(save_dir) / f"Class_distribution_{run_.sample}.png", dpi=600, bbox_inches='tight')

def plot_length_distribution(run_, max_, save_dir=None):
    data_ = pd.DataFrame({
        "length": [len(x) for x in run_.unique["seq_clear"]],
        "seq": run_.unique["seq_clear"]
    })
    
    data_["length"] = data_["length"].astype('category')

    plt.figure(figsize=(4, 3))
    ax = sns.countplot(data=data_, x='length', color='black', width=0.6)

    counts = data_['length'].value_counts().reindex(data_['length'].cat.categories)
    percentages = (counts / counts.sum()) * 100
    annotations = combine(counts.values, percentages.values)
    
    for p, annotation in zip(ax.patches, annotations):
        ax.annotate(annotation, (p.get_x() + p.get_width() / 2., p.get_height()),
                    ha='center', va='bottom', fontsize=8)

    ax.set_ylim(0, max_)
    ax.set_xlabel("Peptide length (aa)")
    ax.set_ylabel("Count")
    ax.set_title(f"{run_.sample} {run_.search_engine} peptides lengths distribution FDR = {run_.fdr}")
    
    ax.yaxis.grid(True, linestyle='--', linewidth=0.5, color='lightgray')
    
    sns.despine()

    if save_dir:
        plt.savefig(Path(save_dir) / f"length_distrib_{run_.sample}_{run_.search_engine}_{run_.fdr}.png", dpi=600, bbox_inches='tight')
    
    for klass in ["Canonical", "Ambiguous", "Non canonical", "Contaminant"]:
        max_for_klass = len(run_.unique[run_.unique['class']==klass])
        data_t = pd.DataFrame({
		"length": [len(x) for x in run_.unique[run_.unique['class']==klass]["seq_clear"]],
		"seq": run_.unique[run_.unique['class']==klass]["seq_clear"]
        })

        data_t["length"] = data_t["length"].astype('category')
        plt.figure(figsize=(4, 3))
        ax = sns.countplot(data=data_t, x='length', color='black', width=0.6)

        counts = data_t['length'].value_counts().reindex(data_t['length'].cat.categories)
        percentages = (counts / counts.sum()) * 100
        annotations = combine(counts.values, percentages.values)

        for p, annotation in zip(ax.patches, annotations):
                ax.annotate(annotation, (p.get_x() + p.get_width() / 2., p.get_height()),
                                        ha='center', va='bottom', fontsize=8)

        ax.set_ylim(0, max_for_klass+50)
        ax.set_xlabel("Peptide length (aa)")
        ax.set_ylabel("Count")
        ax.set_title(f"Total {klass} peptides lengths distribution FDR = {run_.fdr}")

        ax.yaxis.grid(True, linestyle='--', linewidth=0.5, color='lightgray')

        sns.despine()

        if save_dir:
            plt.savefig(Path(save_dir) / f"length_distrib_{klass}_{run_.sample}_{run_.search_engine}_{run_.fdr}.png", dpi=600, bbox_inches='tight')

#rt vs hydrophobycity plotting 
#Source: https://www.biob.in/2014/05/hydrophobicity-plot-using-biopython.html

kd = { 'A': 1.8,'R':-4.5,'N':-3.5,'D':-3.5,'C': 2.5,
       'Q':-3.5,'E':-3.5,'G':-0.4,'H':-3.2,'I': 4.5,
       'L': 3.8,'K':-3.9,'M': 1.9,'F': 2.8,'P':-1.6,
       'S':-0.8,'T':-0.7,'W':-0.9,'Y':-1.3,'V': 4.2 }

palette_kde = {"Canonical": "#a6a6a6", "Non canonical": "#eaaa6a"}

def plot_corr_HI_RT(run_, scale_=kd, save_dir=None):
    HI = []
    rts = []

    data_ = run_.unique[(run_.unique['class']=="Canonical")|(run_.unique['class']=="Non canonical")]

    for index, row in data_.iterrows():
        seq = row['seq_clear']
        num_residue = len(seq)
        values = []

        for residue in seq:
            values.append(scale_[residue])

        HI.append(sum(values))
        rts.append(row['rt'] / 60)

    klass = np.array(data_["class"])
    dataset = pd.DataFrame({'RT': np.array(rts), 'HI': np.array(HI), 'Class': klass}) 

    plt.figure(figsize=(5, 3.5))

    subset_canonical = dataset[dataset['Class'] == "Canonical"]
    sns.scatterplot(data=subset_canonical, x='RT', y='HI', color="#bcbcbc", alpha=0.6, label="Canonical")

    subset_noncanonical = dataset[dataset['Class'] == "Non canonical"]
    sns.scatterplot(data=subset_noncanonical, x='RT', y='HI', color="#eaaa6a", alpha=0.6, label="Non canonical")

    ax = plt.gca()
    for label in dataset['Class'].unique():
        subset = dataset[dataset['Class'] == label]
        sns.kdeplot(x=subset['RT'], y=subset['HI'], ax=ax, fill=False, color=palette_kde[label], levels=5)

    ax.set_xlabel("Retention time (minutes)")
    ax.set_ylabel("Hydrophobicity index")
    ax.set_title(f"{run_.sample} {run_.search_engine} peptides HI vs rt FDR = {run_.fdr}")
    ax.legend()

    ax.yaxis.grid(True, linestyle='--', linewidth=0.5, color='lightgray')
    ax.set(xlim=(0, None))
    sns.despine()

    if save_dir:
        plt.savefig(Path(save_dir) / f"Correlation_HI_RT_{run_.sample}_{run_.search_engine}_{run_.fdr}.png", dpi=600, bbox_inches='tight')

def plot_NetMHCpan_binding_results(run_, xls, save_dir=None):
    with open(xls) as my_file:
        hlas = my_file.readline().rstrip('\n').strip('\t').split('\t\t\t\t')
    df = pd.read_csv(xls, sep="\t", header=1)
    ranks = [column for column in df.columns if column.startswith('EL_Rank')]
    new_cols = ['Peptide', 'NB']
    new_cols.extend(ranks)
    df = df[new_cols].rename(columns=dict(zip(ranks,list(map(lambda x: x+'_rank', hlas)))))
    for hla in hlas:
        df[hla] = df[f'{hla}_rank'].apply(lambda x: 'Strong' if x<0.5 else 'Weak' if (x<2 and x>=0.5) else 'No binder')

    #plot all per allele
    allele_binder_counts = pd.DataFrame.from_dict({'Allele': hlas, 'Count_strong': [len(df[df[hla]=='Strong']) for hla in hlas], 'Count_weak': [len(df[df[hla]=='Weak'])  for hla in hlas]})
    #print(allele_binder_counts)
    bar = allele_binder_counts.plot(kind = 'bar',x='Allele', stacked = True, color = ["#00008b", "#4169e1"], zorder=2)
    bar.yaxis.grid(True, linestyle='--', linewidth=0.5, color='lightgray', zorder=0)
    for i, row in allele_binder_counts.iterrows():
        plt.text(i, row['Count_strong'], str(row['Count_strong']), ha="center", va="bottom", zorder=3)
        plt.text(i, row['Count_strong']+row['Count_weak'], str(row['Count_weak']), ha="center", va="bottom", zorder=3)
    bar.set_ylabel("Count of binders per allele")
    bar.set_title(f"All {run_.sample} peptides \n NetMHCpan binding prediction", fontsize=16, weight='bold')
    plt.xticks(rotation=45)
    sns.despine()
    plt.tight_layout()
    if save_dir:
        plt.savefig(Path(save_dir) /f"{run_.sample}_NetMHCpan_binding_all_per_allele.png", dpi=600, bbox_inches='tight')

    #plot all total
    df['Binder'] = df['NB'] != 0
    plt.figure(figsize=(2, 3))
    bar_plot = sns.barplot(df, y='Binder', errorbar=None, color="#00008b", zorder=3)
    bar_plot.text(0,  df['Binder'].mean(), df['Binder'].sum(), ha="center", va="bottom", zorder=2)
    bar_plot.yaxis.grid(True, linestyle='--', linewidth=0.5, color='lightgray')
    bar_plot.set(ylim=(None, 1.0))
    bar_plot.set_ylabel("Fraction of binders")
    bar_plot.set_title(f"All {run_.sample} peptides \n NetMHCpan binding prediction", fontsize=16, weight='bold')
    sns.despine()
    plt.tight_layout()
    if save_dir:
        plt.savefig(Path(save_dir) /f"{run_.sample}_NetMHCpan_binding_all.png", dpi=600, bbox_inches='tight')

    non_can = run_.unique[run_.unique['class']=='Non canonical']['seq_clear'].tolist()
    mask = df['Peptide'].isin(non_can)
    df_non_can = df[mask]
    #print(df_non_can)
    #plot non canonical per allele
    allele_binder_counts = pd.DataFrame.from_dict({'Allele': hlas, 'Count_strong': [len(df_non_can[df_non_can[hla]=='Strong'])  for hla in hlas], 'Count_weak': [len(df_non_can[df_non_can[hla]=='Weak'])  for hla in hlas]})
    #print(allele_binder_counts)
    bar = allele_binder_counts.plot(kind = 'bar',x='Allele', stacked = True, color = ["#00008b", "#4169e1"], zorder=2)
    bar.yaxis.grid(True, linestyle='--', linewidth=0.5, color='lightgray', zorder=0)
    for i, row in allele_binder_counts.iterrows():
        if row['Count_strong'] != 0:
            plt.text(i, row['Count_strong'], str(row['Count_strong']), ha="center", va="bottom", zorder=3)
        if row['Count_weak'] != 0:
            plt.text(i, row['Count_strong']+row['Count_weak'], str(row['Count_weak']), ha="center", va="bottom", zorder=3)
    bar.set_ylabel("Count of binders per allele")
    bar.set_title(f"Non-canonical peptides \n {run_.sample} \n NetMHCpan binding prediction", fontsize=16, weight='bold')
    plt.xticks(rotation=45)
    sns.despine()
    plt.tight_layout()
    if save_dir:
        plt.savefig(Path(save_dir) /f"{run_.sample}_NetMHCpan_binding_noncanonical_peptides_per_allele.png", dpi=600, bbox_inches='tight')

    #plot non canonical total
    df_non_can['Binder'] = df_non_can['NB'] != 0
    plt.figure(figsize=(2, 3))
    bar_plot = sns.barplot(df_non_can , y='Binder', errorbar=None, color="#00008b", zorder=3)
    bar_plot.text(0,  df_non_can['Binder'].mean(), df_non_can['Binder'].sum(), ha="center", va="bottom", zorder=2)
    bar_plot.yaxis.grid(True, linestyle='--', linewidth=0.5, color='lightgray')
    bar_plot.set(ylim=(None, 1.0))
    bar_plot.set_ylabel("Fraction of binders")
    bar_plot.set_title(f"Non-canonical peptides \n {run_.sample} \n NetMHCpan binding prediction", fontsize=16, weight='bold')
    sns.despine()
    plt.tight_layout()
    if save_dir:
        plt.savefig(Path(save_dir) /f"{run_.sample}_NetMHCpan_binding_noncanonical_peptides_all.png", dpi=600, bbox_inches='tight')

    return df

def main():
    parser = argparse.ArgumentParser(description='Perform postprocessing on Fragpipe outputs')
    parser.add_argument('--input_tsv', required=True, help='Path to Fragpipe output psm.tsv file')
    parser.add_argument('--netMHCpan_xls', required=True, help='Path to NetMHCpan output file')
    parser.add_argument('--GibbsCluster_csv', required=True, help='Path to GibbsCluster output file')
    parser.add_argument('--sample', required=True, help='Project name')
    parser.add_argument('--hla', required=True, help='Sample hlas. Format example: HLA-A02:01,HLA-A03:01,HLA-B27:02,HLA-B44:05,HLA-C02:02')
    parser.add_argument('--output_dir', required=True, help='A directory to store output figures')
    parser.add_argument('--search_engine', required=True, help='MS peptidomics search engine name')
    parser.add_argument('--fdr', required=True, help='MS peptidomics search engine FDR')
    parser.add_argument('--translation_mode', required=True, help='Wildcard translation_mode for file names')
    parser.add_argument('--hla_la_cryptic', required=True, help='A list of cryptic peptides found in benign tissues reported in doi:10.1136/jitc-2020-002071')
    args = parser.parse_args()
    start_time = time.time()
    metadata = pd.DataFrame({'sample': args.sample, 'search_engine': args.search_engine, 'fdr': args.fdr, 'path': args.input_tsv, 'DB': 'uniprot and tumor specific transcripts', 'hla_alleles': args.hla}, index=[0])
    run = MSrun(metadata.iloc[0])
    run.isCanonical()
    psm_count = len(run.psm)
    unique_count = len(run.unique)
    metadata["psm_count"] = psm_count
    metadata["unique_count"] = unique_count
    row_dict = metadata.to_dict(orient="records")[0]
    plot_class_distribution(run, max_ = unique_count+3000, save_dir = args.output_dir)
    plot_length_distribution(run, max_ = 0.6*unique_count+3000,  save_dir = args.output_dir)
    plot_corr_HI_RT(run, scale_=kd, save_dir=args.output_dir)
    df_binding = plot_NetMHCpan_binding_results(run, args.netMHCpan_xls, save_dir = args.output_dir)
    print(df_binding)
    df_final = run.unique.merge(df_binding.rename(columns={'Peptide': 'seq_clear'}), on='seq_clear', how='left')
    columns_to_fillna = dict(zip(args.hla.split(','), ['No binder']*len(args.hla.split(','))))
    columns_to_fillna.update({'Binder': False})
    df_final = df_final.fillna(value=columns_to_fillna)
    hlala_cryptic = pd.read_csv(args.hla_la_cryptic)
    df_final['found_in_benign_tissue'] = df_final[df_final['class']=='Non canonical']['seq_clear'].apply(lambda x: x in hlala_cryptic['Sequence'].values)
    df_final.to_csv(Path(args.output_dir) /f'Summary.{args.sample}.{args.translation_mode}.csv', index=False)

if __name__ == "__main__":
    main()