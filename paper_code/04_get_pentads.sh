for DATASET in schwarzer_TAM schwarzer_NIPBL
do
	python ../pentads/src/get_pentad_cis.py ${DATASET}.cool ./tests_biological/schwarzer_PC1.bed::PC1 --out_pref ./tests_biological/${DATASET}
	python ../pentads/src/plot_pentad.py ./tests_biological/${DATASET}.json --title ${DATASET} --vmin 0.5 --vmax 2 --out_pref ./tests_biological/${DATASET}
	python ../pentads/src/plot_pentad.py ./tests_biological/${DATASET}.json --title ${DATASET} --vmin 0.5 --vmax 2 --format pdf --out_pref ./tests_biological/${DATASET}
	python ../pentads/src/get_pentad_strength.py ${DATASET}.cool ./tests_biological/schwarzer_PC1.bed::PC1 --out_pref ./tests_biological/${DATASET}
done

for DATASET in haarhuis_WT haarhuis_WAPL
do
       python ../pentads/src/get_pentad_cis.py ${DATASET}.cool ./tests_biological/haarhuis_PC1.bed::PC1 --out_pref ./tests_biological/${DATASET}
       python ../pentads/src/plot_pentad.py ./tests_biological/${DATASET}.json --title ${DATASET} --vmin 0.5 --vmax 2 --out_pref ./tests_biological/${DATASET}
       python ../pentads/src/plot_pentad.py ./tests_biological/${DATASET}.json --title ${DATASET} --vmin 0.5 --vmax 2 --format pdf --out_pref ./tests_biological/${DATASET}
	python ../pentads/src/get_pentad_strength.py ${DATASET}.cool ./tests_biological/haarhuis_PC1.bed::PC1 --out_pref ./tests_biological/${DATASET}
done

for DATASET in abramo_0h abramo_2h abramo_3h abramo_4h abramo_5h abramo_8h
do
	python ../pentads/src/get_pentad_cis.py ${DATASET}.cool ./tests_biological/abramo_PC1.bed::PC1 --out_pref ./tests_biological/${DATASET}
	python ../pentads/src/plot_pentad.py ./tests_biological/${DATASET}.json --title ${DATASET} --vmin 0.5 --vmax 2 --out_pref ./tests_biological/${DATASET}
	python ../pentads/src/plot_pentad.py ./tests_biological/${DATASET}.json --title ${DATASET} --vmin 0.5 --vmax 2 --format pdf --out_pref ./tests_biological/${DATASET}
	python ../pentads/src/get_pentad_strength.py ${DATASET}.cool ./tests_biological/abramo_PC1.bed::PC1 --out_pref ./tests_biological/${DATASET}

	python ../pentads/src/get_pentad_distance.py ${DATASET}.cool ./tests_biological/abramo_PC1.bed::PC1 --out_pref ./tests_biological/${DATASET}_distance --distance 10 25 75
	python ../pentads/src/plot_pentad.py ./tests_biological/${DATASET}_distance.json --title ${DATASET}_distance --vmin 0.5 --vmax 2 --closed --out_pref ./tests_biological/${DATASET}_distance
	python ../pentads/src/plot_pentad.py ./tests_biological/${DATASET}_distance.json --title ${DATASET}_distance --vmin 0.5 --vmax 2 --closed --format pdf --out_pref ./tests_biological/${DATASET}_distance
	python ../pentads/src/get_strength_distance.py ${DATASET}.cool ./tests_biological/abramo_PC1.bed::PC1 --out_pref ./tests_biological/${DATASET}_distance --distance 10 25 75
done

for DATASET in du_MII_oocyte du_sperm du_PN5 du_PN5_maternal du_PN5_paternal du_early_2cell du_early_2cell_maternal du_early_2cell_paternal du_late_2cell du_late_2cell_maternal du_late_2cell_paternal du_8cell du_8cell_maternal du_8cell_paternal du_ICM du_ICM_maternal du_ICM_paternal
do
	python ../pentads/src/get_pentad_cis.py ${DATASET}.cool ./tests_biological/du_PC1.bed::PC1 --out_pref ./tests_biological/${DATASET}
	python ../pentads/src/plot_pentad.py ./tests_biological/${DATASET}.json --title ${DATASET} --vmin 0.5 --vmax 2 --out_pref ./tests_biological/${DATASET}
	python ../pentads/src/plot_pentad.py ./tests_biological/${DATASET}.json --title ${DATASET} --vmin 0.5 --vmax 2 --format pdf --out_pref ./tests_biological/${DATASET}
	python ../pentads/src/get_pentad_strength.py ${DATASET}.cool ./tests_biological/du_PC1.bed::PC1 --out_pref ./tests_biological/${DATASET}
done
