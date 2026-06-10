'''
SEARCHES.append(AntigeneSearch(project="example", 
                              input_fasta="path/to/your/neotranscript/fasta", 
                              hla_txt="path/to/your/hla/alleles.txt", 
                              immunopeptidomics_data="path/to/your/immunopeptidomics/data"))

SEARCHES.append(AntigeneSearch(project="SKCM_TAA_MM15",
                                input_fasta="/storage/TAA/_data/estimates/testis/filtered_TAA_testis.fasta",
                                hla_txt="/storage/TAA/_data/IP_skcm/20141208_QEp7_MiBa_SA_HLA-I_MM15_HLA_type.txt",
                                immunopeptidomics_data="/storage/TAA/_data/IP_skcm/20141208_QEp7_MiBa_SA_HLA-I_MM15",
                                translation_mode='stop_to_stop'))
'''
    
from glob import glob 
import os

SEARCHES.append(AntigeneSearch(project="SKCM_TAA_MM15_test_pipeline",
                                input_fasta="/storage/TAA/_data/estimates/testis/filtered_TAA_testis.fasta",
                                hla_txt="/storage/TAA/_data/IPs/IP_skcm/20141208_QEp7_MiBa_SA_HLA-I_MM15_HLA_type.txt",
                                immunopeptidomics_data="/storage/TAA/_data/IPs/IP_skcm/20141208_QEp7_MiBa_SA_HLA-I_MM15",
                                translation_mode='stop_to_stop',
                                grouped_fdr=False))
# paths = [
#     p for p in glob('/potager/_MARIIA/HLA_ligand_atlas_bypatient/*')
#         if os.path.isdir(p) and 'hla_type' not in p] 
# print(paths)

# for path in paths:
#     sample_id = path.split('/')[-1]
#     print(sample_id)
#     SEARCHES.append(
#         AntigeneSearch(
#             project=sample_id.replace('-', '_'),
#             input_fasta="/storage/TAA/_data/estimates/testis/filtered_TAA_testis.fasta",
#             hla_txt=f"/potager/_MARIIA/HLA_ligand_atlas_bypatient/hla_type/{sample_id}_HLA-I.txt",
#             immunopeptidomics_data=path,
#             translation_mode='stop_to_stop'
#         )
#     )