import argparse
import os
import json
import cv2
import numpy as np
import pandas as pd
import multiprocess as mp
import cooler
from cooltools import expected_trans
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
        signal = pd.read_csv(bedgraph_file, header = 0, sep = '\t', names = ['chrom', 'start', 'end', 'PC1'])
        return(list(signal[signal['chrom'] == chromosome]['PC1'].values))
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

def get_area_type(interval_1, interval_2,
                  intervals_A1, intervals_B1,
                  intervals_A2, intervals_B2):
    """
    Find what type of area was extracted from the Hi-C map.

    Parameters:
    interval_1 -- First interval used for extraction.
    interval_2 -- Second interval used for extraction.
    intervals_A1 -- List of compartment A intervals from first chromosome.
    intervals_B1 -- List of compartment B intervals from first chromosome.
    intervals_A2 -- List of compartment A intervals from second chromosome.
    intervals_B2 -- List of compartment B intervals from second chromosome.

    Output:
    Area type.
    """

    if (interval_1 in intervals_A1 and interval_2 in intervals_B2) or\
       (interval_1 in intervals_B1 and interval_2 in intervals_A2):
        return('AB')
    elif (interval_1 in intervals_A1 and interval_2 in intervals_A2):
        return('A')
    elif (interval_1 in intervals_B1 and interval_2 in intervals_B2):
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
parser.add_argument('--max_zeros', default = 0.1, type = float, required = False,
                    help = 'Maximum fraction of bins with zero contacts in an area')
parser.add_argument('--incl_chrms', default='', type = str, required = False,
                    help = 'Chromosomes to include for analysis')
parser.add_argument('--excl_chrms', default = 'Y,M,MT', type = str, required = False,
                    help = 'Chromosomes to exclude from analysis')
parser.add_argument('--out_pref', default = 'pentad_trans', type = str, required = False,
                    help='Prefix for the output files')
parser.add_argument('--save_submatrices', action='store_true', required = False,
                    help='Whether to save extracted submatrices')

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
incl_chrms = args.incl_chrms.split(',')
if incl_chrms != '':
    incl_chrms = incl_chrms + ['chr' + chrm for chrm in incl_chrms]
else:
    pass
excl_chrms = args.excl_chrms.split(',')
excl_chrms = excl_chrms + ['chr' + chrm for chrm in excl_chrms]
out_pref = args.out_pref
save_submatrices = args.save_submatrices

print('Running trans pentad calculation for {} with {}'.format(cool_file, comp_signal))

#############################################
# Read files
#############################################

if not os.path.isfile(cool_file_path):
    raise FileExistsError("cool file with Hi-C matrix doesn't exist")
if not os.path.isfile(comp_signal):
    raise FileExistsError("bedGraph file with compartment signal doesn't exist")

c = cooler.Cooler(cool_file)
chromosomes = c.chroms()[:]['name'].values
if incl_chrms != '':
    chromosomes = [chrm for chrm in chromosomes if chrm in incl_chrms]
chromosomes = [chrm for chrm in chromosomes if chrm not in excl_chrms]
chromsizes = c.chroms()[:][c.chroms()[:]['name'].map(lambda x: x in chromosomes)].set_index('name')
resolution = c.info['bin-size']

#############################################
# Prepare trans expected values
#############################################

trans_df = expected_trans(c)
print('Trans-expected values calculated!')

#############################################
# Calculate average compartment
#############################################

print('Processing trans data...')
average_compartment = [ [], [], [] ]
areas_stats = [ [0], [0], [0] ]

eigenvectors = {}
intervals_A = {}
intervals_B = {}
intervals_zero = {}
all_intervals = {}
for chrm in chromosomes:
    eigenvectors[chrm] = open_eigenvector(comp_signal, chrm, column)
    comp_A_index, comp_B_index, zero_bins = get_compartment_bins(eigenvectors[chrm])
    intervals_A[chrm], intervals_B[chrm], intervals_zero[chrm] = get_compartment_intervals(comp_A_index,
                                                                                           comp_B_index,
                                                                                           zero_bins)
    all_intervals[chrm] = np.sort(intervals_A[chrm] + intervals_B[chrm] + intervals_zero[chrm])

for first_idx in range(len(chromosomes)):

    print('\tChromosome {} trans pairs... '.format(chromosomes[first_idx]))

    for second_idx in range(first_idx+1, len(chromosomes)):

        matrix = c.matrix(balance = True, sparse = True).fetch(chromosomes[first_idx], chromosomes[second_idx]).toarray()
        matrix = np.nan_to_num(matrix) / trans_df[(trans_df['region1'] == chromosomes[first_idx]) & (trans_df['region2'] == chromosomes[second_idx])]['balanced.avg'].values[0]

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

if save_submatrices:
    np.save(out_pref + '.npy', np.array(average_compartment))

average_compartment = [np.nanmedian(x, axis = 0) for x in average_compartment]
print('Average compartment in trans calculated!')

#############################################
# Save output
#############################################

output = {
    'data' : {},
    'stats' : {},
    'type' : 'trans'
}

subplot_titles = ['A', 'B', 'AB']

for area, stat, title in zip(average_compartment, areas_stats, subplot_titles):
    output['data'][title] = area.tolist()
    output['stats'][title] = int(np.sum(stat))

with open(out_pref + '.json', 'w') as w:
    json.dump(output, w)

print('Output saved!')
