'''
SEARCHES.append(AntigeneSearch(project="example", 
                              input_fasta="path/to/your/neotranscript/fasta", 
                              hla_txt="path/to/your/hla/alleles.txt", 
                              immunopeptidomics_data="path/to/your/immunopeptidomics/data"))

SEARCHES.append(AntigeneSearch(project="EPS_DDA_MHCI_GRU1",
                                input_fasta="/mnt/raid0/EPS/EPS_extra_genes.index_specific.groupby_index.fasta",
                                hla_txt="/mnt/raid0/EPS/data/HLAs/GRU_HLAs.txt",
                                immunopeptidomics_data="/mnt/raid0/EPS/data/Immunopeptidomics/DDA_MHC1_GRU1"))

SEARCHES.append(AntigeneSearch(project="example",
                                input_fasta="example/example_MLANA.fasta",
                                hla_txt="example/example_OD5P_hlas.txt",
                                immunopeptidomics_data="_data/example"))                                
'''
SEARCHES.append(AntigeneSearch(project="EPS_DDA_MHCI_GRU1_test_aa",
                                input_aa_fasta="/mnt/raid0/ewing_sarcoma/data/databases/new_db_from_Luuk/20250916_EwS_ncORF.fasta",
                                hla_txt="/mnt/raid0/EPS/data/HLAs/GRU_HLAs.txt",
                                immunopeptidomics_data="/mnt/raid0/EPS/data/Immunopeptidomics/DDA_MHC1_GRU1"))

