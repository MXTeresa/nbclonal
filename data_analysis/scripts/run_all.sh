# this data file has the unfiltered mutations from McCrone et al. eLife 2018 doi.org/10.7554/eLife.35962
# https://github.com/elifesciences-publications/Host_level_IAV_evolution/blob/master/data/processed/secondary/no_cut_trans_freq.csv
iav_dat="data/no_cut_trans_freq.csv"

# this is the metadata file for Popa et al. STM 2020 10.1126/scitranslmed.abe2555
# generated as part of Martin & Koelle STM 2021 10.1126/scitranslmed.abh1803
# https://github.com/koellelab/sarscov2_nb_reanalysis https://github.com/koellelab/sarscov2_nb_reanalysis/blob/main/data_analysis/data/abe255_Data_file_format.csv
sars2_dat='data/abe255_Data_file_format.csv'
# raw data has been re-downloaded and processed according to Martin & Koelle STM 2022 
sars2_vcf='data/seq/*/*_filter_norm.vcf'

# now actually process the iav and sars cov2 data
python3 scripts/process_iav_dat.py \
	--dat $iav_dat

python3 scripts/process_sarscov2_dat.py \
	--metadata $sars2_dat \
	--vcfDir $sars2_vcf 

