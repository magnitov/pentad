# Average compartments calculation for sample data
# Cis
python ./src/get_pentad_cis.py ./test_input/Control.500kb.cool ./test_input/PC1_values_500kb.txt::Control --out_pref ./test_output/Control_cis
python ./src/get_pentad_cis.py ./test_input/Hex5.500kb.cool ./test_input/PC1_values_500kb.txt::Hex5 --out_pref ./test_output/Hex5_cis
# Trans
python ./src/get_pentad_trans.py ./test_input/Control.500kb.cool ./test_input/PC1_values_500kb.txt::Control --out_pref ./test_output/Control_trans
python ./src/get_pentad_trans.py ./test_input/Hex5.500kb.cool ./test_input/PC1_values_500kb.txt::Hex5 --out_pref ./test_output/Hex5_trans
# By distance
python ./src/get_pentad_distance.py ./test_input/Control.500kb.cool ./test_input/PC1_values_500kb.txt::Control --out_pref ./test_output/Control_distance --distance 10 25
python ./src/get_pentad_distance.py ./test_input/Hex5.500kb.cool ./test_input/PC1_values_500kb.txt::Hex5 --out_pref ./test_output/Hex5_distance --distance 10 25

# Visualization of calculated average compartments for sample data
# Cis
python ./src/plot_pentad.py ./test_output/Control_cis.json --title Control --vmin 0.5 --vmax 2 --out_pref ./test_output/Control_cis_plot
python ./src/plot_pentad.py ./test_output/Hex5_cis.json --title Hex5 --vmin 0.5 --vmax 2 --out_pref ./test_output/Hex5_cis_plot
# Trans
python ./src/plot_pentad.py ./test_output/Control_trans.json --title Control --vmin 0.4 --vmax 2.5 --out_pref ./test_output/Control_trans_plot
python ./src/plot_pentad.py ./test_output/Hex5_trans.json --title Hex5 --vmin 0.4 --vmax 2.5 --out_pref ./test_output/Hex5_trans_plot
# By distance
python ./src/plot_pentad.py ./test_output/Control_distance.json --title Control --vmin 0.5 --vmax 2 --out_pref ./test_output/Control_distance_plot
python ./src/plot_pentad.py ./test_output/Hex5_distance.json --title Hex5 --vmin 0.5 --vmax 2 --out_pref ./test_output/Hex5_distance_plot
