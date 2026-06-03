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

SEARCHES.append(AntigeneSearch(project="EPS_DDA_MHCI_GRU1_test_aa",
                                input_aa_fasta="/mnt/raid0/ewing_sarcoma/data/databases/new_db_from_Luuk/20250916_EwS_ncORF.fasta",
                                hla_txt="/mnt/raid0/EPS/data/HLAs/GRU_HLAs.txt",
                                immunopeptidomics_data="/mnt/raid0/EPS/data/Immunopeptidomics/DDA_MHC1_GRU1"))                               

SEARCHES.append(AntigeneSearch(project="EW1260",
                                input_aa_fasta="/mnt/raid0/ewing_sarcoma/data/databases/new_db_from_Luuk/20250916_EwS_ncORF.fasta",
                                data_type="DIA",
                                hla_txt="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1260/EW1260_hlas.txt",
                                immunopeptidomics_data="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1260/mzml"))   

SEARCHES.append(AntigeneSearch(project="EW1343",
                                input_aa_fasta="/mnt/raid0/ewing_sarcoma/data/databases/new_db_from_Luuk/20250916_EwS_ncORF.fasta",
                                data_type="DIA",
                                hla_txt="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1343/EW1343_hlas.txt",
                                immunopeptidomics_data="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1343/mzml"))

SEARCHES.append(AntigeneSearch(project="EW1505",
                                input_aa_fasta="/mnt/raid0/ewing_sarcoma/data/databases/new_db_from_Luuk/20250916_EwS_ncORF.fasta",
                                data_type="DIA",
                                hla_txt="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1505/EW1505_hlas.txt",
                                immunopeptidomics_data="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1505/mzml"))    

SEARCHES.append(AntigeneSearch(project="EW1528",
                                input_aa_fasta="/mnt/raid0/ewing_sarcoma/data/databases/new_db_from_Luuk/20250916_EwS_ncORF.fasta",
                                data_type="DIA",
                                hla_txt="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1528/EW1528_hlas.txt",
                                immunopeptidomics_data="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1528/mzml"))     

SEARCHES.append(AntigeneSearch(project="EW1732",
                                input_aa_fasta="/mnt/raid0/ewing_sarcoma/data/databases/new_db_from_Luuk/20250916_EwS_ncORF.fasta",
                                data_type="DIA",
                                hla_txt="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1732/EW1732_hlas.txt",
                                immunopeptidomics_data="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1732/mzml"))

SEARCHES.append(AntigeneSearch(project="EW1343",
                                input_aa_fasta="/mnt/raid0/ewing_sarcoma/data/databases/new_db_from_Luuk/20250916_EwS_ncORF.fasta",
                                data_type="DIA",
                                hla_txt="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1343/EW1343_hlas.txt",
                                immunopeptidomics_data="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1343/mzml"))

SEARCHES.append(AntigeneSearch(project="EW1505",
                                input_aa_fasta="/mnt/raid0/ewing_sarcoma/data/databases/new_db_from_Luuk/20250916_EwS_ncORF.fasta",
                                data_type="DIA",
                                hla_txt="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1505/EW1505_hlas.txt",
                                immunopeptidomics_data="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1505/mzml"))    
#done
'''
SEARCHES.append(AntigeneSearch(project="EW1528_fdb",
                                input_aa_fasta="/mnt/raid0/ewing_sarcoma/data/databases/new_db_from_Luuk/20260527_EwS_ncORF_filtered.fasta",
                                data_type="DIA",
                                hla_txt="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1528/EW1528_hlas.txt",
                                immunopeptidomics_data="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1528/mzml"))     

SEARCHES.append(AntigeneSearch(project="EW1732_fdb",
                                input_aa_fasta="/mnt/raid0/ewing_sarcoma/data/databases/new_db_from_Luuk/20260527_EwS_ncORF_filtered.fasta",
                                data_type="DIA",
                                hla_txt="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1732/EW1732_hlas.txt",
                                immunopeptidomics_data="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1732/mzml"))

SEARCHES.append(AntigeneSearch(project="EW1260_fdb",
                                input_aa_fasta="/mnt/raid0/ewing_sarcoma/data/databases/new_db_from_Luuk/20260527_EwS_ncORF_filtered.fasta",
                                data_type="DIA",
                                hla_txt="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1260/EW1260_hlas.txt",
                                immunopeptidomics_data="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1260/mzml"))   

SEARCHES.append(AntigeneSearch(project="EW1505_fdb",
                                input_aa_fasta="/mnt/raid0/ewing_sarcoma/data/databases/new_db_from_Luuk/20260527_EwS_ncORF_filtered.fasta",
                                data_type="DIA",
                                hla_txt="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1505/EW1505_hlas.txt",
                                immunopeptidomics_data="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1505/mzml"))    

SEARCHES.append(AntigeneSearch(project="EW1343_fdb",
                                input_aa_fasta="/mnt/raid0/ewing_sarcoma/data/databases/new_db_from_Luuk/20260527_EwS_ncORF_filtered.fasta",
                                data_type="DIA",
                                hla_txt="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1343/EW1343_hlas.txt",
                                immunopeptidomics_data="/mnt/raid0/ewing_sarcoma/data/tumors_immunopeptidomics/EW1343/mzml"))    
