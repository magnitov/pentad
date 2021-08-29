# Pentad

## Introduction

Understanding the role of various factors in 3D genome organisation is essential to determine their impact on shaping large-scale chromatin units such as euchromatin and heterochromatin compartments. However, detailed analysis of changes within and between these compartments is complicated because of a lack of suitable computational methods. To fill this gap, we developed Pentad, a tool to perform calculation, visualisation and quantitative analysis of the average chromatin compartment at different genomic scales.

The average compartment visualisation provided by Pentad represents short- and long-range contacts within A and B compartments together with intercompartmental interactions. The visualisation comprises several types of areas from the Hi-C matrix that are determined based on the annotated A/B compartment signal. Here is the diagram illustrating how the average compartment visualization is created:

![Pentads diagram](https://github.com/magnitov/pentads/blob/development/diagram1.png)

Additionally, Pentad is able calculate compartment strength directly from the average compartments. Here is the diagram demonstrating how this calculations are performed:

![Compartment strength](https://github.com/magnitov/pentads/blob/development/diagram2.png)

## Installation

1. Clone the repo to your local machine:

```
git clone https://github.com/magnitov/pentads.git
cd ./pentads/
```

2. Create conda environment with all the dependencies:

```
conda env create -n pentads-env -f environment.yml
conda activate pentads-env
```

3. Run scripts on the sample dataset (this should take about 15 minutes to complete):

```
bash run_test.sh
```

## Usage

### Parameters

Here is the full description of parameters you may find in our scripts:

*  **cool_file**

Path to the cool file with Hi-C matrix for which you want to create an average compartment.

*  **comp_signal**

Path to the bedGraph file with compartment signal. You could use the ‘::’ syntax to specify a column name from the bedGraph file.

* **rescale_size**

Size to which all the areas in the average compartment will be rescaled. The default is 33 bins.

* **min_dimension**

Minimum dimension of an area (in genomic bins) to be considered for analysis. This is used to filter the short unreliable compartments with noisy compartment signal. The default is 3 bins.

* **max_zeros**

Maximum fraction of bins with zero contacts in an area. This is used to filter sparse regions in the Hi-C matrix. The default is 0.5.

* **cutoff**

Maximum distance between two intervals in the chromosome. We find that compartment regions that are located very far from one another don't add up any significant information, therefore we set a cutoff for the distance between them. This is used only for *cis* and *cis* interactions stratified by distance.

* **distances**

Distance boundaries in Mb separated by space. For example, 10 100 will give <10 Mb, 10-100 Mb, >100 Mb. This is used only for *cis* interactions stratified by distance.

* **excl_chrms**

Chromosomes to exclude from the analysis. By default, we exclude Y,M, and MT, which are sometimes presented in the Hi-C matrices.

* **out_pref**

Prefix for the output files. By default, we save the output to the same directory with prefix *pentad*.

* **center_width**

Fraction of the central fragment of the square that represents intercompartental interactions used for compartment strenght quantification.

* **closed**

For *cis* interactions stratified by distance whether to plot a closed intervals (omitting the last section).

### Average compartment calculation

For calculating the average compartment in any mode you will need two inputs: cool file with Hi-C contact matrix and bedGraph file with compartment signal. The other parameters are optional.

**Average compartment in *cis*:**

```
usage: get_pentad_cis.py [-h] [--rescale_size RESCALE_SIZE]
                         [--min_dimension MIN_DIMENSION]
                         [--max_zeros MAX_ZEROS] [--cutoff CUTOFF]
                         [--excl_chrms EXCL_CHRMS] [--out_pref OUT_PREF]
                         cool_file comp_signal
```

**Average compartment in *trans*:**

```
usage: get_pentad_trans.py [-h] [--rescale_size RESCALE_SIZE]
                           [--min_dimension MIN_DIMENSION]
                           [--max_zeros MAX_ZEROS] [--excl_chrms EXCL_CHRMS]
                           [--out_pref OUT_PREF]
                           cool_file comp_signal
```

**Average compartment in *cis* by distance:**

```
usage: get_pentad_distance.py [-h] [--rescale_size RESCALE_SIZE]
                              [--min_dimension MIN_DIMENSION]
                              [--max_zeros MAX_ZEROS] [--cutoff CUTOFF]
                              --distances DISTANCES [DISTANCES ...]
                              [--excl_chrms EXCL_CHRMS] [--out_pref OUT_PREF]
                              cool_file comp_signal
```

### Visualization

We have a separate script for the visualization of the output files. It automatically detects which type of average compartment was calculated and creates a visualization for you. It requires only the file with calculated pentad, however there are a few parameters to tune your visualization.

```
usage: plot_pentad.py [-h] [--vmin VMIN] [--vmax VMAX] [--cmap CMAP]
                      [--title TITLE] [--closed] [--out_pref OUT_PREF]
                      [--format FORMAT] [--compare COMPARE]
                      average_compartment_path
```

### Compartment strength quantification

For calculating the average compartment strength in any mode you will need same parameters and input files as for average compartment calculation.

**Compartment strength in *cis*:**

```
usage: quant_strength_cis.py [-h] [--rescale_size RESCALE_SIZE]
                             [--min_dimension MIN_DIMENSION]
                             [--max_zeros MAX_ZEROS] [--cutoff CUTOFF]
                             [--center_width CENTER_WIDTH]
                             [--excl_chrms EXCL_CHRMS] [--out_pref OUT_PREF]
                             cool_file comp_signal
```

**Compartment strength in *trans*:**

```
usage: quant_strength_trans.py [-h] [--rescale_size RESCALE_SIZE]
                               [--min_dimension MIN_DIMENSION]
                               [--max_zeros MAX_ZEROS]
                               [--center_width CENTER_WIDTH]
                               [--excl_chrms EXCL_CHRMS] [--out_pref OUT_PREF]
                               cool_file comp_signal
```

**Compartment strength in *cis* by distance:**

```
usage: quant_strength_distance.py [-h] [--rescale_size RESCALE_SIZE]
                                  [--min_dimension MIN_DIMENSION]
                                  [--max_zeros MAX_ZEROS] [--cutoff CUTOFF]
                                  --distances DISTANCES [DISTANCES ...]
                                  [--center_width CENTER_WIDTH]
                                  [--excl_chrms EXCL_CHRMS]
                                  [--out_pref OUT_PREF]
                                  cool_file comp_signal
```

## Citing

If you find this visualization useful for your research, please cite our paper:

Sergey V. Ulianov, Artem K. Velichko, Mikhail D. Magnitov, Artem V. Luzhin, Arkadiy K. Golov, Natalia Ovsyannikova, Igor I. Kireev, Alexey S. Gavrikov, Alexander S. Mishin, Azat K. Garaev, Alexander V. Tyakht, Alexey A. Gavrilov, Omar L. Kantidze, Sergey V. Razin. **Suppression of liquid-liquid phase separation by 1,6-hexanediol partially compromises the 3D genome organization in living cells.** *Nucleic Acids Research* (2021)
* doi: [10.1093/nar/gkab249](https://doi.org/10.1093/nar/gkab249)

## Support
In case you have any questions about the data or need some advice, please contact us via email: mikhail.magnitov@phystech.edu (Mikhail Magnitov) or azatgk@yandex.ru (Azat Garaev).
