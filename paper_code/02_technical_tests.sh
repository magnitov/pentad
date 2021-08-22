# Resolution
for RESOLUTION in 20000 50000 100000 200000 500000 1000000
do
       cooltools call-compartments -o ./tests_technical/rao_GM12878.${RESOLUTION} rao_GM12878.mcool::/resolutions/${RESOLUTION}
done

for RESOLUTION in 1000000 500000 200000 100000 50000 20000
do
	for NUM in {1..10}
	do
		/usr/bin/time -o ./tests_technical/cis_resolution_${RESOLUTION}_test_${NUM}.txt python ../pentads/src/get_pentad_cis.py rao_GM12878.mcool::resolutions/${RESOLUTION} ./tests_technical/rao_GM12878.${RESOLUTION}.cis.vecs.tsv::E1 --out_pref ./tests_technical/out_test_cis
		/usr/bin/time -o ./tests_technical/trans_resolution_${RESOLUTION}_test_${NUM}.txt python ../pentads/src/get_pentad_trans.py rao_GM12878.mcool::resolutions/${RESOLUTION} ./tests_technical/rao_GM12878.${RESOLUTION}.cis.vecs.tsv::E1 --out_pref ./tests_technical/out_test_trans
	done
done

rm ./tests_technical/out_test_cis* ./tests_technical/out_test_trans*

# Compartments size
awk '{ if (NR%2==1) print $0"\t1"; else print $0"\t-1" }' ./tests_technical/hg19.100kb.windows.bed | sed '1s/^/chrom\tstart\tend\tE1\n/' > ./tests_technical/compartments_complexity_100000.tsv
awk '{ if (NR%4==1 || NR%4==2) print $0"\t1"; else print $0"\t-1" }' ./tests_technical/hg19.100kb.windows.bed | sed '1s/^/chrom\tstart\tend\tE1\n/' > ./tests_technical/compartments_complexity_200000.tsv
awk '{ if (NR%10==1 || NR%10==2 || NR%10==3 || NR%10==4 || NR%10==5) print $0"\t1"; else print $0"\t-1" }' ./tests_technical/hg19.100kb.windows.bed | sed '1s/^/chrom\tstart\tend\tE1\n/' > ./tests_technical/compartments_complexity_500000.tsv
awk '{ if (NR%20==1 || NR%20==2 || NR%20==3 || NR%20==4 || NR%20==5 || NR%20==6 || NR%20==7 || NR%20==8 || NR%20==9 || NR%20==10) print $0"\t1"; else print $0"\t-1" }' ./tests_technical/hg19.100kb.windows.bed | sed '1s/^/chrom\tstart\tend\tE1\n/' > ./tests_technical/compartments_complexity_1000000.tsv
awk '{ if (NR%40==1 || NR%40==2 || NR%40==3 || NR%40==4 || NR%40==5 || NR%40==6 || NR%40==7 || NR%40==8 || NR%40==9 || NR%40==10 || NR%40==11 || NR%40==12 || NR%40==13 || NR%40==14 || NR%40==15 || NR%40==16 || NR%40==17 || NR%40==18 || NR%40==19 || NR%40==20) print $0"\t1"; else print $0"\t-1" }' ./tests_technical/hg19.100kb.windows.bed | sed '1s/^/chrom\tstart\tend\tE1\n/' > ./tests_technical/compartments_complexity_2000000.tsv
awk '{ if (NR%100==1 || NR%100==2 || NR%100==3 || NR%100==4 || NR%100==5 || NR%100==6 || NR%100==7 || NR%100==8 || NR%100==9 || NR%100==10 || NR%100==11 || NR%100==12 || NR%100==13 || NR%100==14 || NR%100==15 || NR%100==16 || NR%100==17 || NR%100==18 || NR%100==19 || NR%100==20 || NR%100==21 || NR%100==22 || NR%100==23 || NR%100==24 || NR%100==25 || NR%100==26 || NR%100==27 || NR%100==28 || NR%100==29 || NR%100==30 || NR%100==31 || NR%100==32 || NR%100==33 || NR%100==34 || NR%100==35 || NR%100==36 || NR%100==37 || NR%100==38 || NR%100==39 || NR%100==40 || NR%100==41 || NR%100==42 || NR%100==43 || NR%100==44 || NR%100==45 || NR%100==46 || NR%100==47 || NR%100==48 || NR%100==49 || NR%100==50) print $0"\t1"; else print $0"\t-1" }' ./tests_technical/hg19.100kb.windows.bed | sed '1s/^/chrom\tstart\tend\tE1\n/' > ./tests_technical/compartments_complexity_5000000.tsv

for COMPLEXITY in 5000000 2000000 1000000 500000 200000 100000
do
	for NUM in {1..10}
	do
		/usr/bin/time -o ./tests_technical/cis_complexity_${COMPLEXITY}_test_${NUM}.txt python ../pentads/src/get_pentad_cis.py rao_GM12878.mcool::resolutions/100000 ./tests_technical/compartments_complexity_${COMPLEXITY}.tsv::E1 --out_pref ./tests_technical/out_test_complexity_cis --min_dimension 1
		/usr/bin/time -o ./tests_technical/trans_complexity_${COMPLEXITY}_test_${NUM}.txt python ../pentads/src/get_pentad_trans.py rao_GM12878.mcool::resolutions/100000 ./tests_technical/compartments_complexity_${COMPLEXITY}.tsv::E1 --out_pref ./tests_technical/out_test_complexity_trans --min_dimension 1
        done
done

rm ./tests_technical/out_test_complexity_cis* ./tests_technical/out_test_complexity_trans*
