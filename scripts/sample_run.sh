DD=../empirical_data
COUNTS=sample_peptide_counts.csv
PMETA=peptide_metadata.csv
SMETA=sample_metadata.csv
./PhiP_seq_fold_analysis.py -counts ${DD}/${COUNTS} -peptide_metadata ${DD}/${PMETA} -sample_metadata ${DD}/${SMETA} -corr_th 0.90
