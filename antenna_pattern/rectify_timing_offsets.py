"""
Methods for rectifying for overflows and timing offsets in GNU Radio samples.

Described in https://www.la1k.no/?p=3152.
"""

from combine_samples_and_angles import read_gnuradio_header_element, load_gnuradio_samples
from frequency_analysis import get_center_frequency_and_samplerate
import os
import numpy as np
from scipy.interpolate import interp1d

def get_timestamps_from_header(gnuradio_hdr_file):
    """
    Parse GNU Radio header file for timestamp tags.

    Parameters
    ----------
    gnuradio_hdr_file: str
        Path to GNU Radio header file

    Returns
    -------
    timestamps: ndarray
        List over timestamps
    num_items: ndarray
        List over number of items for which the timestamps are valid
    """

    handle = open(gnuradio_hdr_file)

    #get number of header elements contained in header file
    _, element_size = read_gnuradio_header_element(handle)
    handle.seek(0, 0)

    header_size = os.path.getsize(gnuradio_hdr_file)
    num_elements = header_size/element_size

    #get all header information
    timestamps = np.zeros(num_elements)
    samplerates = np.zeros(num_elements)
    num_items = np.zeros(num_elements)
    num_read = 0

    for i in np.arange(0, num_elements):
        info, element_size = read_gnuradio_header_element(handle)
        if element_size == 0:
            break

        num_read += element_size
        handle.seek(num_read, 0)

        timestamps[i] = info['rx_time']
        num_items[i] = info['nitems']

    return timestamps, num_items

def copy_first_part_of_header(filename_in, filename_out):
    """
    Copy first part of header file to a new header file.

    Parameters
    ----------
    filename_in: str
        Input filename. Input header is assumed to be filename_in.hdr.
    filename_out: str
        Output filename. Output header is assumed to be filename_out.hdr.
    """
    header_file_out = open(filename_out + '.hdr', 'w')
    header_file_in = open(filename_in + '.hdr')
    _, header_size = read_gnuradio_header_element(header_file_in)
    header_file_in.seek(0)
    header_data = header_file_in.read(header_size)

    header_file_out.write(header_data)

def fill_overflows(filename_in, filename_out):
    """
    Takes in a GNU Radio meta file sink file and fills places with overflows
    with a number of samples corresponding to the missing samples during the
    overflow in order to rectify somewhat for the timing issues.

    Parameters
    ----------
    filename_in: str
        Input filename
    filename_out: str
        Output filename for writing the filled samples. Will also produce a
        header file, containing the first header element of the original file.
    """

    timestamps, num_items = get_timestamps_from_header(filename_in + '.hdr')
    f_c, f_s = get_center_frequency_and_samplerate(filename_in)
    _, data = load_gnuradio_samples(filename_in, return_full_timestamps=False)

    #calculate missing samples
    expected_samples = np.diff(timestamps)*f_s
    missing_samples = expected_samples - num_items[0:-1]
    missing_samples = np.round(missing_samples).astype(int)

    #find indices where there is overflow problems
    problem_samples = np.arange(0, len(missing_samples))[missing_samples > 0]
    sample_counter = np.cumsum(num_items).astype(int)

    #fill missing samples, write to new file
    out_file = open(filename_out, 'w')
    prev_problem = 0
    for problem_sample in problem_samples:
        #contiguous part of the data
        out_file.write(data[sample_counter[prev_problem]:sample_counter[problem_sample]].tobytes())

        #add a number of samples corresponding to the number of missing samples, interpolated between the current and next sample
        interpolated_samples = np.repeat(interp1d([0, 1], [data[problem_sample], data[problem_sample+1]])(0.5), missing_samples[problem_sample]).astype(np.complex64)

        out_file.write(interpolated_samples.tobytes())

        prev_problem = problem_sample

    #write out last part of the data
    out_file.write(data[sample_counter[problem_sample]:].tobytes())

    #write header
    copy_first_part_of_header(filename_in, filename_out)

import sys
if __name__ == "__main__":
    if len(sys.argv) <= 2:
        print("Rectifies for overflows in GNU Radio meta file sink files by filling in the missing samples.\n")
        print("Usage: " + sys.argv[0] + " input_filename output_filename")
    else:
        fill_overflows(sys.argv[1], sys.argv[2])
