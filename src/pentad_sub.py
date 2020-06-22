import argparse
import os
import h5py
import cv2
import numpy as np
import pandas as pd
import cooler
from cooltools import numutils
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

def get_subcompartment_bins(eigenvector, subs):
    compartment = {}
    for sub in subs:
        compartment[sub] = [ind for (ind, eig) in zip(np.arange(len(eigenvector)), eigenvector) if eig == sub]
    return(compartment)

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

def get_area_type(interval_1, interval_2, comp_intervals, subs):
    type_of_area = []
    for sub1 in subs:
        if interval_1 in comp_intervals[sub1]:
            type_of_area.append(sub1)
            for sub2 in subs:
                if interval_2 in comp_intervals[sub2]:
                    type_of_area.append(sub2)

    type_of_area.sort()
    if type_of_area[0] == type_of_area[1]:
        return(type_of_area[0])
    else:
        return(f'{type_of_area[0]}-{type_of_area[1]}')

def area_dimensions_are_large_enough(img, min_dimension):
    return(len(img) >= min_dimension and len(img[0]) >= min_dimension)

def area_has_enough_data(img, max_zeros_fraction):
    return(len(np.where(img.ravel() == 0)[0]) < max_zeros_fraction * len(img.ravel()))

def area_is_close_enough(intervals_1, intervals_2, matrix_size, cutoff):
    return(np.mean(intervals_2) < np.mean(intervals_1) + cutoff * matrix_size)

def area_is_at_diagonal(i, j):
    return(i == j)


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
parser.add_argument('--excl_chrms', default='Y,M,MT', type = str, required = False,
                    help = 'Chromosomes to exclude from analysis')
# Plot
parser.add_argument('--vmin', default = 0.5, type = float, required = False,
                    help = 'Lower limit for the colormap')
parser.add_argument('--vmax', default = 2, type = float, required = False,
                    help = 'Upper limit for the colormap')
parser.add_argument('--cmap', default = 'coolwarm', type = str, required = False,
                    help = 'Colormap to use for the visualization')
parser.add_argument('--title', default = '', type = str, required = False,
                    help = 'Suptitle to use for the visualization')
# Output
parser.add_argument('--out_pref', default = 'pentad', type = str, required = False,
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
resolution = c.info['bin-size']
subs = ['A1', 'A2', 'B1', 'B2', 'B3']

# Make list of subplot titles
subplot_titles = [f'{sub} short' for sub in subs]
subplot_titles += [f'{sub} long' for sub in subs]
for i in range(0, 5):
    for j in range(i+1,5):
        subplot_titles.append(f'{subs[i]}-{subs[j]}')

# Calculate average compartment
average_compartment = { i : [] for i in subplot_titles}
areas_stats = { i : [0] for i in subplot_titles}
for chromosome in chromosomes:
    print('Chromosome {}...'.format(chromosome))

    eigenvector = open_eigenvector(comp_signal, chromosome)
    comp_index = get_subcompartment_bins(eigenvector, subs)
    comp_intervals = { sub : calculate_intervals_from_range(comp_index[sub]) for sub in subs}
    all_intervals = np.concatenate(tuple([comp_intervals[sub] for sub in subs]))
    all_intervals.sort()

    matrix = c.matrix(balance = True, sparse = True).fetch(chromosome).toarray()
    matrix = np.nan_to_num(matrix)
    matrix, *other = numutils.observed_over_expected(matrix)

    for i in range(0, len(all_intervals)):
        for j in range(i, len(all_intervals)):
            area = get_area_from_matrix(matrix, all_intervals[i], all_intervals[j])

            if area_dimensions_are_large_enough(area, min_dimension) and\
                area_has_enough_data(area, max_zeros) and\
                area_is_close_enough(all_intervals[i], all_intervals[j],
                                     len(matrix), distance_cutoff):

                area_resized = resize_area(area, rescale_size)
                area_type = get_area_type(all_intervals[i], all_intervals[j],
                                          intervals_A, intervals_B)

                if len(area_type) == 2:
                    if area_is_at_diagonal(i, j):
                        average_compartment[f'{area_type} short'].append(area_resized)
                    else:
                        average_compartment[f'{area_type} long'].append(area_resized)
                else:
                    average_compartment[area_type].append(area_resized)


    for i in subplot_titles:
        areas_stats[i].append(len(average_compartment[i])-np.sum(areas_stats[i]))


average_compartment = { x: np.nanmedian(average_compartment[x], axis = 0) for x in average_compartment.keys()}

print('Average compartment calculated!')
print('Total areas calculated:')
for i in range(5):
    print(f'\t{subplot_titles[i]}: {np.sum(areas_stats[subplot_titles[i]])}\
          \t{subplot_titles[i+5]}: {np.sum(areas_stats[subplot_titles[i+5]])}')
for sbplt in subplot_titles[10:]:
    print(f'\t{sbplt}: {np.sum(areas_stats[sbplt])}')

# Visualize average compartment
fig = plt.figure(figsize = ( 20, 24 ))
plt.suptitle(title, x = 0.5125, y = 0.925, fontsize = 22)
for i in range(5):
    plt.subplot(6,5,i+1)
    plt.imshow(average_compartment[subplot_titles[i]], cmap = cmap, norm = LogNorm(vmax = vmax, vmin = vmin))
    plt.title(subplot_titles[i], fontsize = 20)
    plt.xticks([], [])
    plt.yticks([], [])

    plt.subplot(6,5,i+6)
    plt.imshow(average_compartment[subplot_titles[i+5]], cmap = cmap, norm = LogNorm(vmax = vmax, vmin = vmin))
    plt.title(subplot_titles[i+5], fontsize = 20)
    plt.xticks([], [])
    plt.yticks([], [])

indices = [12, 13, 14, 15, 18, 19, 20, 24, 25, 30]
for i in range(10):
    plt.subplot(6,5,indices[i])
    plt.imshow(average_compartment[subplot_titles[i+10]], cmap = cmap, norm = LogNorm(vmax = vmax, vmin = vmin))
    plt.title(subplot_titles[i+10], fontsize = 20)
    plt.xticks([], [])
    plt.yticks([], [])

cbar_ax = fig.add_axes([0.95, 0.25, 0.02, 0.5])
cbar = plt.colorbar(cax = cbar_ax)

plt.savefig(out_pref + '.png', bbox_inches = 'tight')
plt.clf()

print('Visualization created!')
