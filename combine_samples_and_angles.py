"""
Functions for loading angles produced using rotctld_angle_printer.py, and
reading GNU Radio samples produced using a file metasink with detached header,
and for interpolating angles based on the GNU Radio timestamps.

Usage example:

```
import combine_samples_and_angles as comb
angles = comb.read_angles('angle_file')
timestamps, samples = comb.load_gnuradio_samples('gnuradio_file')

timestamps, power_spectrum = comb.power_spectrum(timestamps, samples)

azimuth = comb.interpolate_azimuth(timestamps, angles)

#if frequency bin number 236 measures the frequency of the signal we want to estimate the power of:
plt.plot(azimuth, power_spectrum[:, 236])
```
"""

import pandas as pd
import numpy as np
import scipy.interpolate
import gnuradio.blocks as blocks
from gnuradio.blocks import parse_file_metadata
import gnuradio.gr as gr

def load_angles(angle_file):
    """
    Read angles and timestamps from file
    generated using rotctld_angle_printer.py.

    \\param angle_file File containing angles and corresponding timestamps
    \\return Dataframe containing angles and timestamps
    """
    angles = pd.read_table(angle_file)
    angles.timestamp = pd.to_datetime(angles.timestamp)
    return angles

def spectrum(timestamps, iq_samples, start_bin=0, num_bins=None, n_fft=8192, return_timestamps=True):
    """
    Calculate FFT spectrum of a given signal.

    \\param timestamps, for changing the rate of the timestamps to the rate of
    the FFT windows
    \\param iq_samples Raw IQ samples as obtained from GNU Radio
    \\param start_bin Start bin for the bins we want to keep
    \\param num_bins Number of bins we want to keep. If None, will return all bins
    \\param n_fft Size of FFT
    \\param return_timestamps Whether to return resampled timestamps
    \\return Resampled timestamps (if to be returned), FFT spectrum
    """

    if num_bins is None:
        num_bins = n_fft

    num_frames = len(iq_samples)/n_fft
    end_bin = start_bin + num_bins

    ret_spectrum = np.zeros((num_frames, num_bins), dtype=np.complex64)
    for i in np.arange(0, num_frames):
        fft_res = np.fft.fftshift(np.fft.fft(iq_samples[i*n_fft:(i+1)*n_fft]))
        ret_spectrum[i,:] = fft_res[start_bin:end_bin]

    if return_timestamps:
        timestamps = timestamps[n_fft/2:-1:n_fft]
        timestamps = timestamps[0:ret_spectrum.shape[0]]
        return timestamps, ret_spectrum
    else:
        return ret_spectrum

def spectrum_to_db(spectrum):
    return 10*np.log10(np.abs(spectrum)**2)

def power_spectrum(timestamps, iq_samples, start_bin = 0, num_bins = None, n_fft = 8192, return_timestamps=True):
    """
    Calculate the power spectrum of a given signal.

    For the antenna diagrams, we plot a specific frequency bin as a function of
    azimuth angle.

    See spectrum() for arguments.

    \\return Resampled timestamps, power spectrum
    """

    ret_values = spectrum(timestamps, iq_samples, start_bin, num_bins, n_fft, return_timestamps)

    if return_timestamps:
        return ret_values[0], spectrum_to_db(ret_values[1])
    else:
        return spectrum_to_db(ret_values)

import pmt

def load_gnuradio_header(gnuradio_hdr_file):
    """
    Load GNU Radio meta file header. Function load_gnuradio_samples()
    uses this function to read the file header before reading the file data,
    it is not necessary to call it explicitly.

    The header file will probably contain multiple header instances,
    one for each issue of a new tag, but we read only the first.

    \\param gnuradio_hdr_file GNU Radio header filename
    \\return Header info
    """
    handle = open(gnuradio_hdr_file)
    header_str = handle.read(parse_file_metadata.HEADER_LENGTH)
    header = pmt.deserialize_str(header_str)
    info = parse_file_metadata.parse_header(header, False)

    #get extra information
    if info["extra_len"] > 0:
        extra_str = handle.read(info["extra_len"])
        extra = pmt.deserialize_str(extra_str)
        extra_info = parse_file_metadata.parse_extra_dict(extra, info, False)

    return info

def load_gnuradio_samples(gnuradio_file, return_full_timestamps=True):
    """
    Read gnuradio samples and corresponding timestamps from file.

    File should be produced using a file meta sink, with detached_header=True.
    If vectors are output, or the stream is decimated, make sure that
    the relative rate change of the file meta sink is set to the correct
    rate. This function will assume that the rx_rate (sample rate) specifies
    the sample rate of each vector, corresponding to a row in the output
    data matrix.

    It is assumed that the source sends a correct rx_time-tag so that this
    corresponds to a UNIX timestamp (USRP does this).

    \\param gnuradio_file Filename
    \\param return_full_timestamps Whether to construct and return full set of timestamps for each sample (True), or just the timestamp for first and last sample (False)
    \\return Tuple of timestamps and gnuradio data
    samples matrix (float, samples x vec_length).

    If return_full_timestamps is set to true, the timestamps
    will be a vector of length samples x 1. Otherwise, timestamps
    is of length 2 x 1 and will contain the first and last timestamp.
    """

    #read file header
    header = load_gnuradio_header(gnuradio_file + ".hdr")
    if header['cplx']:
        datatype = np.complex64
    else:
        datatype = np.float32
    vec_length = header['size']/datatype(1).itemsize

    #read in data
    data = np.memmap(gnuradio_file, offset=0, dtype=datatype, mode='r')

    if (vec_length > 1):
        data = np.reshape(data, (-1, vec_length))

    first_timestamp = np.datetime64(pd.to_datetime(header['rx_time'], unit='s'))
    sample_period_ns = 1.0/header['rx_rate']*1.0e09

    if return_full_timestamps:
        #construct timestamps, assuming constant sample rate and no dropped samples
        seconds_since_beginning = np.arange(0, np.shape(data)[0])*sample_period_ns * np.timedelta64(1, 'ns')
        timestamp = first_timestamp + seconds_since_beginning
    else:
        timestamp = [first_timestamp, first_timestamp + np.shape(data)[0]*sample_period_ns*np.timedelta64(1, 'ns')]
    return timestamp, data

def interpolate_angles(interpolation_timestamps, angles, direction='azimuth'):
    """
    Interpolate azimuth/elevation angles with respect to timestamps. Find which azimuth/elevation
    angles the interpolation timestamps correspond.

    The first angle timestamp is counted as the "epoch", and the timestamps are
    recalculated as the number of seconds from this. Samples outside the angle
    timerange are counted to belong to the first and last azimuth angle.

    \\param interpolation_timestamps Timestamps on which to interpolate
    \\angles Angles dataframe as read using load_angles(), containing
    timestamps and azimuth angles as columns
    \\direction Whether to yield azimuth ('azimuth') or elevation ('elevation')
    \\return Azimuth or elevation angles in array of length corresponding to the input
    timestamps
    """
    first_timestamp = angles.timestamp.as_matrix()[0].copy()
    last_timestamp = angles.timestamp.as_matrix()[-1].copy()

    subset = (interpolation_timestamps >= first_timestamp) & (interpolation_timestamps <= last_timestamp)
    subset_timestamps = interpolation_timestamps[subset].copy()

    #set timestamps to be number of seconds since first timestamp
    angles_timestamp = (angles.timestamp - first_timestamp).as_matrix()/np.timedelta64(1, 's')
    subset_timestamps = (subset_timestamps - first_timestamp)/np.timedelta64(1, 's')

    #interpolate recorded angles with respect to timestamps of subset samples
    angle_interpolator = scipy.interpolate.interp1d(angles_timestamp, angles[direction].as_matrix())
    interpolated_angle = angle_interpolator(subset_timestamps)

    #add the samples that were taken outside of the timestamp ranges of the angle measurements
    first_part = np.repeat(angles[direction].iloc[0], np.count_nonzero(interpolation_timestamps < first_timestamp))
    last_part = np.repeat(angles[direction].iloc[-1], np.count_nonzero(interpolation_timestamps > last_timestamp))
    return np.concatenate((first_part, interpolated_angle, last_part))
