import argparse
import os
import h5py
import cv2
import numpy as np
import pandas as pd
import multiprocess as mp
import cooler
from cooltools import numutils
from cooltools.expected import trans_expected, blocksum_pairwise
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
matplotlib.use('Agg')
import seaborn as sns
sns.set_context('poster')
import warnings
warnings.filterwarnings('ignore')

# Compartment signal processing
def open_eigenvector(bedgraph_file, chromosome):
    signal = pd.read_csv(bedgraph_file, header = None, sep = '\t')
    return(list(signal[signal[0] == chromosome][3].values))

def get_compartment_bins(eigenvector):
    compartment_A = [ind for (ind, eig) in zip(np.arange(len(eigenvector)), eigenvector) if eig > 0]
    compartment_B = [ind for (ind, eig) in zip(np.arange(len(eigenvector)), eigenvector) if eig < 0]
    zero_bins = [ind for (ind, eig) in zip(np.arange(len(eigenvector)), eigenvector) if eig == 0]
    return(compartment_A, compartment_B, zero_bins)

def calculate_intervals_from_range(list_range):
    list_range = list(list_range)
    intervals = []
    for idx, item in enumerate(list_range):
        if not idx or item-1 != intervals[-1][-1]:
            intervals.append([item])
        else:
            intervals[-1].append(item)
    return(intervals)

def get_compartment_intervals(compartment_A, compartment_B, zero_bins):
    intervals_A = calculate_intervals_from_range(compartment_A)
    intervals_B = calculate_intervals_from_range(compartment_B)
    intervals_zero = calculate_intervals_from_range(zero_bins)
    return(intervals_A, intervals_B, intervals_zero)

# Hi-C map areas processing
def get_area_from_matrix(matrix, intervals_list_1, intervals_list_2):
    return(matrix[np.ix_(intervals_list_1, intervals_list_2)])

def resize_area(img, bin_size):
    img_resized = cv2.resize(img * 255 / max(img.ravel()), (bin_size, bin_size))
    img_resized = img_resized / 255 * max(img.ravel())
    return(img_resized)

def get_area_type(interval_1, interval_2,
                  intervals_A1, intervals_B1,
                  intervals_A2, intervals_B2):
    if (interval_1 in intervals_A1 and interval_2 in intervals_B2) or\
       (interval_1 in intervals_B1 and interval_2 in intervals_A2):
        return('AB')
    elif (interval_1 in intervals_A1 and interval_2 in intervals_A2):
        return('A')
    elif (interval_1 in intervals_B1 and interval_2 in intervals_B2):
        return('B')

def area_dimensions_are_large_enough(img, min_dimension):
    return(len(img) >= min_dimension and len(img[0]) >= min_dimension)

def area_has_enough_data(img, max_zeros_fraction):
    return(len(np.where(img.ravel() == 0)[0]) < max_zeros_fraction * len(img.ravel()))


# PARSING ARGUMENTS
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Main arguments
parser.add_argument('cool_file', type = str,
                    help = 'Path to the cool file with Hi-C matrix')
parser.add_argument('comp_signal', type = str,
                    help = 'Path to the bedGraph file with compartment signal')
# Extra control parameters
parser.add_argument('--rescale_size', default = 33, type = int, required = False,
                    help = 'Size to rescale all areas in average compartment')
parser.add_argument('--min_dimension', default = 3, type = int, required = False,
                    help = 'Minimum dimension of an area (in genomic bins)')
parser.add_argument('--max_zeros', default = 0.1, type = float, required = False,
                    help = 'Maximum fraction of bins with zero contacts in an area')
parser.add_argument('--distance', default = 0.75, type = float, required = False,
                    help = 'Maximum distance between two intervals in chromosome fractions')
parser.add_argument('--excl_chrms', default = 'Y,M,MT', type = str, required = False,
                    help = 'Chromosomes to exclude from analysis')

# Plot
parser.add_argument('--vmin', default = 0.4, type = float, required = False,
                    help = 'Lower limit for the colormap')
parser.add_argument('--vmax', default = 2.5, type = float, required = False,
                    help = 'Upper limit for the colormap')
parser.add_argument('--cmap', default = 'coolwarm', type = str, required = False,
                    help = 'Colormap to use for the visualization')
parser.add_argument('--title', default = '', type = str, required = False,
                    help = 'Suptitle to use for the visualization')
# Output
parser.add_argument('--out_pref', default = 'pentad_trans', type = str, required = False,
                    help='Prefix for the output files')

args = parser.parse_args()


# Parse arguments
cool_file = args.cool_file
comp_signal = args.comp_signal

rescale_size = args.rescale_size
min_dimension = args.min_dimension
max_zeros = args.max_zeros
distance_cutoff = args.distance
excl_chrms = args.excl_chrms.split(',')
excl_chrms = excl_chrms + ['chr' + chrm for chrm in excl_chrms]

vmin = args.vmin
vmax = args.vmax
cmap = args.cmap
title = args.title

out_pref = args.out_pref

# Check that cool and bedGraph files exist
if not os.path.isfile(cool_file):
    raise FileExistsError("cool file with Hi-C matrix doesn't exist")
if not os.path.isfile(comp_signal):
    raise FileExistsError("bedGraph file with compartment signal doesn't exist")

# Read Hi-C matrix
c = cooler.Cooler(cool_file)
chromosomes = c.chroms()[:]['name'].values
chromosomes = [chrm for chrm in chromosomes if chrm not in excl_chrms]
chromsizes = c.chroms()[:][c.chroms()[:]['name'].map(lambda x: x in chromosomes)].set_index('name')
resolution = c.info['bin-size']

# Calculate expected values
supports = [(chrm, 0, chromsizes.loc[chrm,'length']) for chrm in chromosomes]
balanced_transform = {"balanced": lambda pixels: pixels["count"] * pixels["weight1"] * pixels["weight2"]}

records = blocksum_pairwise(c, supports,
                            transforms=balanced_transform,
                            chunksize=resolution,
                            map=mp.Pool().map)

trans_records = {
        ( region1[0], region2[0] ): val for ( region1, region2 ), val in records.items()
        }

trans_df = pd.DataFrame.from_dict( trans_records, orient="index" )
trans_df.index.rename( ["chrom1","chrom2"], inplace=True )
trans_df["balanced.avg"] = trans_df["balanced.sum"] / trans_df["n_valid"]
print('Trans-expected values calculated')

# Calculate average compartment
average_compartment = [ [], [], [] ]
areas_stats = [ [0], [0], [0] ]

eigenvectors = {}
intervals_A = {}
intervals_B = {}
intervals_zero = {}
all_intervals = {}
for chrm in chromosomes:
    eigenvectors[chrm] = open_eigenvector(comp_signal, chrm)
    comp_A_index, comp_B_index, zero_bins = get_compartment_bins(eigenvectors[chrm])
    intervals_A[chrm], intervals_B[chrm], intervals_zero[chrm] = get_compartment_intervals(comp_A_index,
                                                                                           comp_B_index,
                                                                                           zero_bins)
    all_intervals[chrm] = np.sort(intervals_A[chrm] + intervals_B[chrm] + intervals_zero[chrm])

for first_idx in range(len(chromosomes)):

    print('Chromosome {}... '.format(chromosomes[first_idx]))

    for second_idx in range(first_idx+1, len(chromosomes)):

        matrix = c.matrix(balance = True, sparse = True).fetch(chromosomes[first_idx], chromosomes[second_idx]).toarray()
        matrix = np.nan_to_num(matrix) / trans_df.loc[chromosomes[first_idx], chromosomes[second_idx]]['balanced.avg']

        for i in range(len(all_intervals[chromosomes[first_idx]])):
            for j in range(len(all_intervals[chromosomes[second_idx]])):
                area = get_area_from_matrix(matrix,
                                            all_intervals[chromosomes[first_idx]][i],
                                            all_intervals[chromosomes[second_idx]][j])

                if area_dimensions_are_large_enough(area, min_dimension) and\
                area_has_enough_data(area, max_zeros):

                    area_resized = resize_area(area, rescale_size)
                    area_type = get_area_type(all_intervals[chromosomes[first_idx]][i], all_intervals[chromosomes[second_idx]][j],
                                              intervals_A[chromosomes[first_idx]], intervals_B[chromosomes[first_idx]],
                                              intervals_A[chromosomes[second_idx]], intervals_B[chromosomes[second_idx]])
                    if area_type == 'A':
                        average_compartment[0].append(area_resized)
                    elif area_type == 'B':
                        average_compartment[1].append(area_resized)
                    elif area_type == 'AB':
                        average_compartment[2].append(area_resized)

        for i in range(3):
            areas_stats[i].append(len(average_compartment[i])-np.sum(areas_stats[i]))

average_compartment = [np.nanmedian(x, axis = 0) for x in average_compartment]
np.save(out_pref + '.npy', np.array(average_compartment))
print('Average compartment calculated!')
print('Total areas piled-up:\n\tA: {}\n\tB: {}\
                            \n\tbetween A and B: {}'.format(np.sum(areas_stats[0]),
                                                            np.sum(areas_stats[1]),
                                                            np.sum(areas_stats[2])))

# Visualize average compartment
subplot_titles = ['A', 'B', 'AB']
subplot_indexes = [1, 2, 3]

fig = plt.figure(figsize = (12, 4))
plt.suptitle(title, x = 0.5125, y = 1.02, fontsize = 22)

for layout, subtitle, index in zip(average_compartment, subplot_titles, subplot_indexes):
    plt.subplot(1, 3, index)
    plt.title(subtitle, fontsize = 15)
    plt.imshow(layout, cmap = cmap, norm = LogNorm(vmax = vmax, vmin = vmin))
    plt.xticks([], [])
    plt.yticks([], [])

cbar_ax = fig.add_axes([0.95, 0.15, 0.02, 0.7])
cbar = plt.colorbar(cax = cbar_ax)

plt.savefig(out_pref + '.png', bbox_inches = 'tight')
plt.clf()

print('Visualization created!')
