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

def get_area_type(interval_1, interval_2, intervals_A, intervals_B):
    if (interval_1 in intervals_A and interval_2 in intervals_B) or\
       (interval_1 in intervals_B and interval_2 in intervals_A):
        return('AB')
    elif (interval_1 in intervals_A and interval_2 in intervals_A):
        return('A')
    elif (interval_1 in intervals_B and interval_2 in intervals_B):
        return('B')

def area_dimensions_are_large_enough(img, min_dimension):
    return(len(img) >= min_dimension and len(img[0]) >= min_dimension)

def area_has_enough_data(img, max_zeros_fraction):
    return(len(np.where(img.ravel() == 0)[0]) < max_zeros_fraction * len(img.ravel()))

def area_is_close_enough(intervals_1, intervals_2, matrix_size, cutoff):
    return(np.mean(intervals_2) < np.mean(intervals_1) + cutoff * matrix_size)


def get_distance_index(intervals_1, intervals_2, distance_intervals, resolution):
    index = 0
    distance = (np.mean(intervals_2) - np.mean(intervals_1)) * resolution
    while distance > distance_intervals[index % len(distance_intervals)] and index < len(distance_intervals):
        index += 1
    return(index)

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
parser.add_argument('--cutoff', default = 0.75, type = float, required = False,
                    help = 'Maximum distance between two intervals in chromosome fractions')
parser.add_argument('--distances', nargs = '+', type = int, required = True,
                    help = 'Distance boundaries for <> in Mb') # Finish!  Maybe use string.split(',') instead of nargs
                    # maybe [0] will work!
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
distance_cutoff = args.cutoff
distance_intervals = [i*10**6 for i in args.distances]
interval_number = len(distance_intervals) + 1
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

# Calculate average compartment
average_compartment = [ [ [], [], [] ] for i in range(interval_number)]
areas_stats = [ [ [0], [0], [0] ] for i in range(interval_number)]
for chromosome in chromosomes:
    print('Chromosome {}...'.format(chromosome))

    eigenvector = open_eigenvector(comp_signal, chromosome)
    comp_A_index, comp_B_index, zero_bins = get_compartment_bins(eigenvector)
    intervals_A, intervals_B, intervals_zero = get_compartment_intervals(comp_A_index,
                                                                         comp_B_index,
                                                                         zero_bins)
    all_intervals = np.sort(intervals_A + intervals_B + intervals_zero)

    matrix = c.matrix(balance = True, sparse = True).fetch(chromosome).toarray()
    matrix = np.nan_to_num(matrix)
    matrix, *other = numutils.observed_over_expected(matrix)

    for i in range(0, len(all_intervals)):
        for j in range(i, len(all_intervals)):
            if not area_is_at_diagonal(i, j):
                area = get_area_from_matrix(matrix, all_intervals[i], all_intervals[j])

                if area_dimensions_are_large_enough(area, min_dimension) and\
                    area_has_enough_data(area, max_zeros) and\
                    area_is_close_enough(all_intervals[i], all_intervals[j],
                                         len(matrix), distance_cutoff):

                    area_resized = resize_area(area, rescale_size)
                    area_type = get_area_type(all_intervals[i], all_intervals[j],
                                          intervals_A, intervals_B)

                    index = get_distance_index(all_intervals[i], all_intervals[j], distance_intervals, resolution)

                    if area_type == 'A':
                        average_compartment[index][0].append(area_resized)
                    elif area_type == 'B':
                        average_compartment[index][1].append(area_resized)
                    elif area_type == 'AB':
                        average_compartment[index][2].append(area_resized)

    for i in range(interval_number):
        for j in range(3):
            areas_stats[i][j].append(len(average_compartment[i][j])-np.sum(areas_stats[i][j]))

for i in range(interval_number):
    average_compartment[i] = [np.nanmedian(x, axis = 0) for x in average_compartment[i]]

np.save(out_pref + '.npy', np.array(average_compartment))
print('Average compartment calculated!')

distance_titles = []
distance_titles.append('>{} Mb'.format(distance_intervals[-1]/10**6))
if interval_number > 1:
    distance_titles.insert(0,'<{} Mb'.format(distance_intervals[0]/10**6))

if interval_number > 2:
    for i in range(interval_number - 2):
        distance_titles.insert(i+1,'{}–{} Mb'.format(distance_intervals[i]/10**6,
                                                     distance_intervals[i+2]/10**6))

print('Total areas piled-up:')
for dist_title in distance_titles:
    print('\t{}:\n\tA: {}\n\t\tB: {}\
                         \n\t\tAB: {}\n'.format(dist_title,
                                              np.sum(areas_stats[0]),
                                              np.sum(areas_stats[1]),
                                              np.sum(areas_stats[2])))


# Visualize average compartment
row_titles = ['A', 'B', 'AB']
subplot_indexes = []
for i in range(interval_number):
    subplot_indexes += [ i + 1,
                         i + 1 + interval_number,
                         i + 1 + 2 * interval_number]

fig = plt.figure(figsize = ( interval_number * 4, 12 ))
#plt.suptitle(title, x = 0.5125, y = 0.98, fontsize = 22)

for layout, index in zip(average_compartment, subplot_indexes):
    plt.subplot(interval_number, 3, index)
    plt.imshow(layout, cmap = cmap, norm = LogNorm(vmax = vmax, vmin = vmin))
    plt.xticks([], [])
    plt.yticks([], [])

cbar_ax = fig.add_axes([0.95, 0.25, 0.02, 0.5])
cbar = plt.colorbar(cax = cbar_ax)

plt.savefig(out_pref + '.png', bbox_inches = 'tight')
plt.clf()

print('Visualization created!')
