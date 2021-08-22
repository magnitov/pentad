JUICER="/home/magnitov/Software/juicer_tools_1.22.01.jar"

# Rao et al. 2017
wget -O rao_GM12878.hic https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_insitu_primary%2Breplicate_combined_30.hic
hic2cool convert -r 5000 rao_GM12878.hic rao_GM12878.5kb.cool
cooler zoomify -p 16 --balance --balance-args "--max-iters 500" -r 5000,10000,20000,50000,100000,200000,500000,1000000,2000000,5000000 -o rao_GM12878.mcool rao_GM12878.5kb.cool
rm rao_GM12878.5kb.cool rao_GM12878.hic

# Schwarzer et al. 2017
wget -O schwarzer_TAM.cool.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93431/suppl/GSE93431_TAM.100kb.cool.HDF5.gz
wget -O schwarzer_NIPBL.cool.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE93nnn/GSE93431/suppl/GSE93431_NIPBL.100kb.cool.HDF5.gz
gunzip schwarzer*.gz

# Haarhuis et al. 2017
wget -O haarhuis_WT.allValidPairs.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE95nnn/GSE95014/suppl/GSE95014_Hap1.validPairs.txt.gz
wget -O haarhuis_WAPL.allValidPairs.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE95nnn/GSE95014/suppl/GSE95014_WaplKO_3.3.validPairs.txt.gz
gunzip haarhuis*.gz

for DATASET in haarhuis_WT haarhuis_WAPL
do
	awk '{ print $4"\t"$2"\t"$3"\t"$7"\t"$5"\t"$6 }' ${DATASET}.allValidPairs.txt | sed 's/+/0/g' | sed 's/-/1/g' | awk '{ if ($2>$5) print $4"\t"$5"\t"$6"\t0\t"$1"\t"$2"\t"$3"\t1"; else print $1"\t"$2"\t"$3"\t0\t"$4"\t"$5"\t"$6"\t1" }' | sort -k2,2d -k6,6d > ${DATASET}.pre
	java -jar ${JUICER} pre -j 8 -r 100000 -k KR ${DATASET}.pre ${DATASET}.hic hg19
	hic2cool convert -r 100000 ${DATASET}.hic ${DATASET}.cool
	cooler balance --max-iters 500 ${DATASET}.cool
	rm ${DATASET}.allValidPairs.txt ${DATASET}.pre ${DATASET}.hic
done

# Abramo et al. 2019
wget -O abramo_0h_rep1.mcool.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3909nnn/GSM3909703/suppl/GSM3909703_TB-HiC-Dpn-R2-T0_hg19.1000.multires.cool.gz
wget -O abramo_2h_rep1.mcool.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3909nnn/GSM3909698/suppl/GSM3909698_TB-HiC-Dpn-R2-T2_hg19.1000.multires.cool.gz
wget -O abramo_3h_rep1.mcool.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3909nnn/GSM3909694/suppl/GSM3909694_TB-HiC-Dpn-R2-T3_hg19.1000.multires.cool.gz
wget -O abramo_4h_rep1.mcool.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3909nnn/GSM3909691/suppl/GSM3909691_TB-HiC-Dpn-R2-T4_hg19.1000.multires.cool.gz
wget -O abramo_5h_rep1.mcool.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3909nnn/GSM3909689/suppl/GSM3909689_TB-HiC-Dpn-R2-T5_hg19.1000.multires.cool.gz
wget -O abramo_8h_rep1.mcool.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3909nnn/GSM3909686/suppl/GSM3909686_TB-HiC-Dpn-R2-T8_hg19.1000.multires.cool.gz
wget -O abramo_0h_rep2.mcool.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3909nnn/GSM3909722/suppl/GSM3909722_TB-HiC-Dpn-T0_hg19.1000.multires.cool.gz
wget -O abramo_2h_rep2.mcool.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3909nnn/GSM3909718/suppl/GSM3909718_TB-HiC-Dpn-T2_hg19.1000.multires.cool.gz
wget -O abramo_3h_rep2.mcool.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3909nnn/GSM3909716/suppl/GSM3909716_TB-HiC-Dpn-T3_hg19.1000.multires.cool.gz
wget -O abramo_4h_rep2.mcool.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3909nnn/GSM3909714/suppl/GSM3909714_TB-HiC-Dpn-T4_hg19.1000.multires.cool.gz
wget -O abramo_5h_rep2.mcool.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3909nnn/GSM3909712/suppl/GSM3909712_TB-HiC-Dpn-T5_hg19.1000.multires.cool.gz
wget -O abramo_8h_rep2.mcool.gz https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3909nnn/GSM3909709/suppl/GSM3909709_TB-HiC-Dpn-T8_hg19.1000.multires.cool.gz
gunzip abramo*.gz

for DATASET in 0h 2h 3h 4h 5h 8h
do
	cooler coarsen -k 100 -o abramo_${DATASET}_rep1.cool abramo_${DATASET}_rep1.mcool::/14
	cooler coarsen -k 100 -o abramo_${DATASET}_rep2.cool abramo_${DATASET}_rep2.mcool::/14
	cooler merge abramo_${DATASET}.cool abramo_${DATASET}_rep1.cool abramo_${DATASET}_rep2.cool
	cooler balance --max-iters 500 abramo_${DATASET}.cool
done
rm abramo*.mcool abramo*rep*.cool

# Du at al. 2017
wget -O du_MII_oocyte.allValidPairs.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE82nnn/GSE82185/suppl/GSE82185_MII_rep12_allValidPairs.txt.gz
wget -O du_sperm.allValidPairs.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE82nnn/GSE82185/suppl/GSE82185_sperm_rep123_allValidPairs.txt.gz
wget -O du_PN5.allValidPairs.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE82nnn/GSE82185/suppl/GSE82185_PN5_rep1234_allValidPairs.txt.gz
wget -O du_PN5_maternal.allValidPairs.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE82nnn/GSE82185/suppl/GSE82185_PN5_rep1234_maternal_allValidPairs.txt.gz
wget -O du_PN5_paternal.allValidPairs.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE82nnn/GSE82185/suppl/GSE82185_PN5_rep1234_paternal_allValidPairs.txt.gz
wget -O du_early_2cell.allValidPairs.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE82nnn/GSE82185/suppl/GSE82185_early_2cell_rep123_allValidPairs.txt.gz
wget -O du_early_2cell_maternal.allValidPairs.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE82nnn/GSE82185/suppl/GSE82185_early_2cell_rep123_maternal_allValidPairs.txt.gz
wget -O du_early_2cell_paternal.allValidPairs.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE82nnn/GSE82185/suppl/GSE82185_early_2cell_rep123_paternal_allValidPairs.txt.gz
wget -O du_late_2cell.allValidPairs.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE82nnn/GSE82185/suppl/GSE82185_late_2cell_rep1234_allValidPairs.txt.gz
wget -O du_late_2cell_maternal.allValidPairs.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE82nnn/GSE82185/suppl/GSE82185_late_2cell_rep1234_maternal_allValidPairs.txt.gz
wget -O du_late_2cell_paternal.allValidPairs.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE82nnn/GSE82185/suppl/GSE82185_late_2cell_rep1234_paternal_allValidPairs.txt.gz
wget -O du_8cell.allValidPairs.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE82nnn/GSE82185/suppl/GSE82185_8cell_rep123_allValidPairs.txt.gz
wget -O du_8cell_maternal.allValidPairs.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE82nnn/GSE82185/suppl/GSE82185_8cell_rep123_maternal_allValidPairs.txt.gz
wget -O du_8cell_paternal.allValidPairs.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE82nnn/GSE82185/suppl/GSE82185_8cell_rep123_paternal_allValidPairs.txt.gz
wget -O du_ICM.allValidPairs.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE82nnn/GSE82185/suppl/GSE82185_ICM_rep123_allValidPairs.txt.gz
wget -O du_ICM_maternal.allValidPairs.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE82nnn/GSE82185/suppl/GSE82185_ICM_rep123_maternal_allValidPairs.txt.gz
wget -O du_ICM_paternal.allValidPairs.txt.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE82nnn/GSE82185/suppl/GSE82185_ICM_rep123_paternal_allValidPairs.txt.gz
gunzip *.gz

for DATASET in du_MII_oocyte du_sperm du_PN5 du_PN5_maternal du_PN5_paternal du_early_2cell du_early_2cell_maternal du_early_2cell_paternal du_late_2cell du_late_2cell_maternal du_late_2cell_paternal du_8cell du_8cell_maternal du_8cell_paternal du_ICM du_ICM_maternal du_ICM_paternal
do
	awk '{ print $4"\t"$2"\t"$3"\t"$7"\t"$5"\t"$6 }' ${DATASET}.allValidPairs.txt | sed 's/+/0/g' | sed 's/-/1/g' | awk '{ if ($2>$5) print $4"\t"$5"\t"$6"\t0\t"$1"\t"$2"\t"$3"\t1"; else print $1"\t"$2"\t"$3"\t0\t"$4"\t"$5"\t"$6"\t1" }' | sort -k2,2d -k6,6d > ${DATASET}.pre
	java -jar ${JUICER} pre -j 8 -r 100000 -k KR ${DATASET}.pre ${DATASET}.hic mm9
	hic2cool convert -r 100000 ${DATASET}.hic ${DATASET}.cool
	cooler balance --max-iters 500 ${DATASET}.cool
	rm ${DATASET}.allValidPairs.txt ${DATASET}.pre ${DATASET}.hic
done
