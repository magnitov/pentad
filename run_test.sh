# Download sample data
wget -O Ulianov_et_al_Control.mcool https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138543/suppl/GSE138543_Control.mcool
wget -O Ulianov_et_al_Hex5.mcool https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138543/suppl/GSE138543_Hex5.mcool
wget -O Ulianov_et_al_Recovery.mcool https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138543/suppl/GSE138543_Recovery.mcool
wget -O compartment_signal.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE138nnn/GSE138543/suppl/GSE138543_compartment_signal.100kb.txt.gz
gunzip compartment_signal.txt.gz

# Average compartments calculation for sample data
for DATASET in Control Hex5 Recovery
do
	# Cis
	python ./src/get_pentad_cis.py Ulianov_et_al_${DATASET}.mcool::resolutions/100000 compartment_signal.txt::${DATASET} --out_pref ./test/${DATASET}_cis
	# Trans
	python ./src/get_pentad_trans.py Ulianov_et_al_${DATASET}.mcool::resolutions/100000 compartment_signal.txt::${DATASET} --max_zeros 0.5 --out_pref ./test/${DATASET}_trans
	# By distance
	python ./src/get_pentad_distance.py Ulianov_et_al_${DATASET}.mcool::resolutions/100000 compartment_signal.txt::${DATASET} --out_pref ./test/${DATASET}_distance --distance 10 25
done

# Visualization of calculated average compartments for sample data
for DATASET in Control Hex5 Recovery
do
	# Cis
	python ./src/plot_pentad.py ./test/${DATASET}_cis.json --title ${DATASET} --out_pref ./test/${DATASET}_cis_plot
	# Trans
	python ./src/plot_pentad.py ./test/${DATASET}_trans.json --title ${DATASET} --vmin 0.05 --vmax 20 --out_pref ./test/${DATASET}_trans_plot
	# By distance
	python ./src/plot_pentad.py ./test/${DATASET}_distance.json --title ${DATASET} --out_pref ./test/${DATASET}_distance_plot
done

# Average compartment strength quantification
for DATASET in Control Hex5 Recovery
do
	# Cis
	python ./src/quant_strength_cis.py Ulianov_et_al_${DATASET}.mcool::resolutions/100000 compartment_signal.txt::${DATASET} --out_pref ./test/${DATASET}_cis_strength
	# Trans
	python ./src/quant_strength_trans.py Ulianov_et_al_${DATASET}.mcool::resolutions/100000 compartment_signal.txt::${DATASET} --max_zeros 0.5 --out_pref ./test/${DATASET}_trans_strength
	# By distance
	python ./src/quant_strength_distance.py Ulianov_et_al_${DATASET}.mcool::resolutions/100000 compartment_signal.txt::${DATASET} --out_pref ./test/${DATASET}_distance_strength --distance 10 25
done

# Remove sample data
rm *.mcool compartment_signal.txt
