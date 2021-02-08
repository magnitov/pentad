import argparse
import os
import json
import cv2
import numpy as np
import pandas as pd
import cooler
from cooltools import numutils
import warnings
warnings.filterwarnings('ignore')

#############################################
# Compartment signal processing functions
#############################################

def open_eigenvector(bedgraph_file, chromosome, column = None):
    """
    Get an eigenvector for specific chromosome from a BED-like file.

    Parameters:
    bedgraph_file -- BED-like file from cooltools call-compartments function with eigenvectors.
    chromosome -- Name of a chromosome for which eigenvector should be extracted.
    column -- Name of a column in bedgraph file that contains eigenvector values.

    Output:
    A list with eigenvector corresponding to specified chromosome.
    """

    if column == None:
        signal = pd.read_csv(bedgraph_file, header = None, sep = '\t')
        return(list(signal[signal['chrom'] == chromosome][3].values))
    else:
        signal = pd.read_csv(bedgraph_file, header = 0, sep = '\t')
        return(list(signal[signal['chrom'] == chromosome][column].values))

#############################################

def get_compartment_bins(eigenvector):
    """
    Get indexes of compartments A/B and NaN (zero) bins from eigenvector.

    Parameters:
    eigenvector -- A list with eigenvector values for specific chromosome.

    Output:
    Three lists with indexes of compartment A/B and zeros bins.
    """

    compartment_A = [ind for (ind, eig) in zip(np.arange(len(eigenvector)), eigenvector) if eig > 0]
    compartment_B = [ind for (ind, eig) in zip(np.arange(len(eigenvector)), eigenvector) if eig < 0]
    zero_bins = [ind for (ind, eig) in zip(np.arange(len(eigenvector)), eigenvector) if eig == 0]
    return(compartment_A, compartment_B, zero_bins)

#############################################

def calculate_intervals_from_range(list_range):
    """
    Merge a list of continuous indexes into list of intervals.

    Parameters:
    list_range -- A list with compartment indexes to be merged into intervals.

    Output:
    A list of lists with intervals.
    """

    list_range = list(list_range)
    intervals = []
    for idx, item in enumerate(list_range):
        if not idx or item-1 != intervals[-1][-1]:
            intervals.append([item])
        else:
            intervals[-1].append(item)
    return(intervals)

#############################################

def get_compartment_intervals(compartment_A, compartment_B, zero_bins):
    """
    Apply intervals merging to compartment A/B and zeros bins.

    Parameters:
    compartment_A -- indexes of eigenvector corresponding to compartment A bins.
    compartment_B -- indexes of eigenvector corresponding to compartment B bins.
    zero_bins -- indexes of eigenvector corresponding to zero bins.

    Output:
    Intervals of compartment A/B and zeros bins.
    """

    intervals_A = calculate_intervals_from_range(compartment_A)
    intervals_B = calculate_intervals_from_range(compartment_B)
    intervals_zero = calculate_intervals_from_range(zero_bins)
    return(intervals_A, intervals_B, intervals_zero)

#############################################
# Hi-C map areas processing functions
#############################################

def get_area_from_matrix(matrix, intervals_list_1, intervals_list_2):
    """
    Extract an area from Hi-C map that lies at the intersection of two intervals.

    Parameters:
    matrix -- Full Hi-C map of a chromosome.
    intervals_list_1 -- First interval for extraction.
    intervals_list_2 -- Second interval for extraction.

    Output:
    Extracted Hi-C map area.
    """

    return(matrix[np.ix_(intervals_list_1, intervals_list_2)])

#############################################

def resize_area(img, bin_size):
    """
    Rescale Hi-C map area to a square with defined edge size.

    Parameters:
    img -- Hi-C map area extracted.
    bin_size -- Rescaling sequare size in bins.

    Output:
    Rescaled Hi-C map area.
    """

    img_resized = cv2.resize(img * 255 / max(img.ravel()), (bin_size, bin_size))
    img_resized = img_resized / 255 * max(img.ravel())
    return(img_resized)

#############################################

def get_area_type(interval_1, interval_2, intervals_A, intervals_B):
    """
    Find what type of area was extracted from the Hi-C map.

    Parameters:
    interval_1 -- First interval used for extraction.
    interval_2 -- Second interval used for extraction.
    intervals_A -- List of compartment A intervals.
    intervals_B -- List of compartment B intervals.

    Output:
    Area type.
    """

    if (interval_1 in intervals_A and interval_2 in intervals_B) or\
       (interval_1 in intervals_B and interval_2 in intervals_A):
        return('AB')
    elif (interval_1 in intervals_A and interval_2 in intervals_A):
        return('A')
    elif (interval_1 in intervals_B and interval_2 in intervals_B):
        return('B')

#############################################

def area_dimensions_are_large_enough(img, min_dimension):
    """
    Find whether extracted area is large enough for the analysis.

    Parameters:
    img -- Hi-C map area extracted.
    min_dimensions -- Minimum size for area dimensions.

    Output:
    Whether an area is large enough.
    """

    return(len(img) >= min_dimension and len(img[0]) >= min_dimension)

#############################################

def area_has_enough_data(img, max_zeros_fraction):
    """
    Find whether extracted area has enough data for the analysis.

    Parameters:
    img -- Hi-C map area extracted.
    max_zeros_fraction -- Maximum fraction of zeros in area.

    Output:
    Whether an area has enough data.
    """

    return(len(np.where(img.ravel() == 0)[0]) < max_zeros_fraction * len(img.ravel()))

#############################################

def area_is_close_enough(intervals_1, intervals_2, matrix_size, cutoff):
    """
    Find whether extracted area has enough data for the analysis.

    Parameters:
    intervals_1 -- First interval used for extraction.
    intervals_2 -- Second interval used for extraction.
    matrix_size -- Size of full Hi-C map of chromosome.
    cutoff -- Maximum distance between intervals as chromosome size fraction.

    Output:
    Whether an area is close enough to the diagonal.
    """

    return(np.mean(intervals_2) < np.mean(intervals_1) + cutoff * matrix_size)

#############################################

def area_is_at_diagonal(i, j):
    """
    Find whether extracted area is from the diagonal.

    Parameters:
    i -- First interval index.
    j -- Second interval index.

    Output:
    Whether an area is at the diagonal.
    """

    return(i == j)

#############################################
# Parse arguments
#############################################

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('cool_file', type = str,
                    help = 'Path to the cool file with Hi-C matrix')
parser.add_argument('comp_signal', type = str,
                    help = 'Path to the bedGraph file with compartment signal. Use the \‘::\’ syntax to specify a column name')
parser.add_argument('--rescale_size', default = 33, type = int, required = False,
                    help = 'Size to rescale all areas in average compartment')
parser.add_argument('--min_dimension', default = 3, type = int, required = False,
                    help = 'Minimum dimension of an area (in genomic bins)')
parser.add_argument('--max_zeros', default = 0.5, type = float, required = False,
                    help = 'Maximum fraction of bins with zero contacts in an area')
parser.add_argument('--cutoff', default = 0.75, type = float, required = False,
                    help = 'Maximum distance between two intervals in chromosome fractions')
parser.add_argument('--excl_chrms', default='Y,M,MT', type = str, required = False,
                    help = 'Chromosomes to exclude from analysis')
parser.add_argument('--out_pref', default = 'pentad', type = str, required = False,
                    help='Prefix for the output files')

args = parser.parse_args()

cool_file = args.cool_file
if '::' in cool_file:
    cool_file_path = cool_file[:cool_file.find('::')]
else:
    cool_file_path = cool_file
comp_signal = args.comp_signal.split('::')
if len(comp_signal) == 2:
    column = comp_signal[1]
else:
    column = None
comp_signal = comp_signal[0]

rescale_size = args.rescale_size
min_dimension = args.min_dimension
max_zeros = args.max_zeros
distance_cutoff = args.cutoff
excl_chrms = args.excl_chrms.split(',')
excl_chrms = excl_chrms + ['chr' + chrm for chrm in excl_chrms]
out_pref = args.out_pref

print('Running cis pentad calculation for {} with {}'.format(cool_file, comp_signal))

#############################################
# Read files
#############################################

if not os.path.isfile(cool_file_path):
    raise FileExistsError("cool file with Hi-C matrix doesn't exist")
if not os.path.isfile(comp_signal):
    raise FileExistsError("bedGraph file with compartment signal doesn't exist")

c = cooler.Cooler(cool_file)
chromosomes = c.chroms()[:]['name'].values
chromosomes = [chrm for chrm in chromosomes if chrm not in excl_chrms]
resolution = c.info['bin-size']

#############################################
# Calculate average compartment
#############################################

print('Processing cis data...')
average_compartment = [[], [], [], [], []]
areas_stats = [[0], [0], [0], [0], [0]]
for chromosome in chromosomes:
    print('\tChromosome {}...'.format(chromosome))

    eigenvector = open_eigenvector(comp_signal, chromosome, column)
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
            area = get_area_from_matrix(matrix, all_intervals[i], all_intervals[j])

            if area_dimensions_are_large_enough(area, min_dimension) and\
                area_has_enough_data(area, max_zeros) and\
                area_is_close_enough(all_intervals[i], all_intervals[j],
                                     len(matrix), distance_cutoff):

                area_resized = resize_area(area, rescale_size)
                area_type = get_area_type(all_intervals[i], all_intervals[j],
                                          intervals_A, intervals_B)

                if area_is_at_diagonal(i, j):
                    if area_type == 'A':
                        average_compartment[0].append(area_resized)
                    elif area_type == 'B':
                        average_compartment[1].append(area_resized)
                else:
                    if area_type == 'A':
                        average_compartment[2].append(area_resized)
                    elif area_type == 'B':
                        average_compartment[3].append(area_resized)
                    elif area_type == 'AB':
                        average_compartment[4].append(area_resized)

    for i in range(0, 5):
        areas_stats[i].append(len(average_compartment[i])-np.sum(areas_stats[i]))


average_compartment = [np.nanmedian(x, axis = 0) for x in average_compartment]
print('Average compartment in cis calculated!')

#############################################
# Save output
#############################################

output = {
    'data' : {},
    'stats' : {},
    'type' : 'cis'
}

subplot_titles = ['Short-range A', 'Short-range B',
                  'Long-range A', 'Long-range B',
                  'Between A and B']

for area, stat, title in zip(average_compartment, areas_stats, subplot_titles):
    output['data'][title] = area.tolist()
    output['stats'][title] = int(np.sum(stat))

with open(out_pref + '.json', 'w') as w:
    json.dump(output, w)

print('Output saved!')
