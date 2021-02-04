# pentads

## Introduction

Understanding the effect of various factors on the 3D genome organization is essential to determine their roles in shaping large-scale chromatin structure units such as euchromatin (A) and heterochromatin (B) compartments. However, currently available tools for the pile-up analysis of the Hi-C maps can only deal with chromatin domains and loops. Therefore, we developed a simple tool to compute and visualize the average Hi-C compartment. 

It is currently implemented in three python scripts for (i) averaging cis interactions, (ii) averaging trans interactions, and (iii) averaging cis interactions with stratification by genomic distance. Here is the diagram illustrating how the average compartment visualization is created:

![Pentads diagram](https://github.com/magnitov/pentads/blob/master/diagram.png)

This visualization (we named it pentad) is aimed to represent short- and long-range contacts within A and B compartments and contacts between A and B compartments. It consists of several piled-up areas from the observed-over-expected Hi-C matrix that are determined based on the compartment signal provided by the user. The tool can filter the areas based on their dimensions in genomic bins, the amount of zero-contact pixels, and the distance between the bases that form the area. The areas that pass the filters are extracted from the matrix and rescaled using bilinear interpolation. Rescaled areas of the same type are averaged genome-wide and aggregated into one pentad. 

## Installation

1. Create conda environment with all the dependencies:

```
conda env create -n pentads-env -f environment.yml 
conda activate pentads-env
```

2. Clone the repo to your local machine:

```
git clone https://github.com/magnitov/pentads.git
```

3. Run scripts on the test dataset (this should take a couple of minutes to complete):

```
cd ./pentads/
bash run_test.sh
```

## Usage

### Calculations

For calculating the average compartment in any mode you will need two inputs: cool file with Hi-C contact matrix and bedGraph file with compartment signal. The other parameters are optional. Here is the full description of parameters you may find in our scripts:

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

Maximum distance between two intervals in the chromosome. We find that compartment regions that are located very far from one another don't add up any significant information, therefore we set a cutoff for the distance between them. This is used only for *cis* interactions.

* **distances**

Distance boundaries in Mb separated by space. For example, 10 100 will give <10 Mb, 10-100 Mb, >100 Mb. This is used only for *cis* interactions stratified by distance.

* **excl_chrms**

Chromosomes to exclude from the analysis. By default, we exclude Y,M, and MT, which are sometimes presented in the Hi-C matrices.

* **out_pref**

Prefix for the output files. By default, we save the output to the same directory with prefix *pentad*.

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
                      [--format FORMAT]
                      average_compartment_path
```

## Citing

If you find this visualization useful for your research, please cite our paper:

Sergey V. Ulianov, Artem K. Velichko, Mikhail D. Magnitov, Artem V. Luzhin, Arkadiy K. Golov, Natalia Ovsyannikova, Igor I. Kireev, Alexey S. Gavrikov, Alexander S. Mishin, Azat K. Garaev, Alexander V. Tyakht, Alexey A. Gavrilov, Omar L. Kantidze, Sergey V. Razin. **Suppression of liquid-liquid phase separation by 1,6-hexanediol partially compromises the 3D genome organization in living cells.** *bioRxiv* (2020)
* doi: [10.1101/2020.05.18.101261](https://doi.org/10.1101/2020.05.18.101261)

## Support
In case you have any questions or need some advice, please contact us via email: mikhail.magnitov@phystech.edu (Mikhail Magnitov) or azatgk@yandex.ru (Azat Garaev).
