# Average compartments calculation for sample data
# Cis
python ./src/get_pentad_cis.py ./test/Control.100kb.cool ./test/PC1_values.txt::Control --out_pref ./out/Control_cis
python ./src/get_pentad_cis.py ./test/Hex5.100kb.cool ./test/PC1_values.txt::Hex5 --out_pref ./out/Hex5_cis
# Trans
python ./src/get_pentad_trans.py ./test/Control.100kb.cool ./test/PC1_values.txt::Control --out_pref ./out/Control_trans
python ./src/get_pentad_trans.py ./test/Hex5.100kb.cool ./test/PC1_values.txt::Hex5 --out_pref ./out/Hex5_trans
# By distance
python ./src/get_pentad_distance.py ./test/Control.100kb.cool ./test/PC1_values.txt::Control --out_pref ./out/Control_distance --distance 10 50 100
python ./src/get_pentad_distance.py ./test/Hex5.100kb.cool ./test/PC1_values.txt::Hex5 --out_pref ./out/Hex5_distance --distance 10 50 100

# Visualization of calculated average compartments for sample data
# Cis
python ./src/plot_pentad.py ./out/Control_cis.json --title Control --out_pref ./out/Control_cis_plot
python ./src/plot_pentad.py ./out/Hex5_cis.json --title Hex5 --out_pref ./out/Hex5_cis_plot
# Trans
python ./src/plot_pentad.py ./out/Control_trans.json --title Control --out_pref ./out/Control_trans_plot
python ./src/plot_pentad.py ./out/Hex5_trans.json --title Hex5 --out_pref ./out/Hex5_trans_plot
# By distance
python ./src/plot_pentad.py ./out/Control_distance.json --title Control --out_pref ./out/Control_distance_plot
python ./src/plot_pentad.py ./out/Hex5_distance.json --title Hex5 --out_pref ./out/Hex5_distance_plot
