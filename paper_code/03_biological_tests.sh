mkdir -p tests_biological

for DATASET in schwarzer_TAM haarhuis_WT du_ICM 
do
	cooltools call-compartments -o ./tests_biological/${DATASET} ${DATASET}.cool
done
