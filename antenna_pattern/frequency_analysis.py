"""
Helper functions for frequency analysis and picking/interpolation of frequency
bins to yield a time-variant signal for a specific frequency of interest.

See interpolate_at_frequency(filename, frequency, n_fft=1024) for an entry method.
"""

from combine_samples_and_angles import load_gnuradio_header, load_gnuradio_samples, load_angles, interpolate_angles
import numpy as np
import pmt
from numpy.fft import fftshift, fft
from scipy.signal import blackman
from scipy.interpolate import interp1d

def interpolate_at_frequency(filename, frequency, n_fft=1024):
    """
    Read samples from file, calculate FFT and interpolate exactly at the specified frequency using cubic splines.

    Parameters
    ----------
    filename : str
        GNU Radio filename
    frequency :
        Frequency of interest
    n_fft : optional
        Length of FFT

    Returns
    -------
    stamps :
        Timestamps for first and last sample
    samples :
        FFT samples at frequency of interest
    """

    delta_f = 5000 #extract frequency bins through frequency \pm delta_f
    freqs, fft, stamps = load_and_calculate_spectrum(filename, frequency=frequency, n_fft=n_fft, delta_f = delta_f)

    samples = interp1d(freqs, np.abs(fft), kind='cubic')(frequency)
    return stamps, samples

def load_and_calculate_spectrum(gnuradio_filename, frequency, delta_f = 0, n_fft=8192, num_seconds=None):
    """
    Load IQ samples from file and calculate the FFT spectrum around a desired frequency.

    Parameters
    ----------
    gnuradio_filename : str
        Path to GNU Radio file
    frequency :
        Center frequency of interest in Hz
    delta_f : optional
        Frequency range of returned spectrum
    n_fft : optional
        Length of FFT
    num_seconds : optional
        Number of seconds to include in the returned spectrum.
    Returns
    -------
    freqs :
        Frequency axis
    fft :
        Discrete fourier transform of the input signal at specified time and frequency ranges (num_timestamps x num_frequencies)
    timestamps :
        Time axis
    """

    #frequencies for each fft bin
    f_c, f_s = get_center_frequency_and_samplerate(gnuradio_filename)
    freqs = frequency_axis(f_c, f_s, n_fft)
    freq_resolution = freqs[1] - freqs[0]

    #get desired frequency subset
    delta_f = np.max([freq_resolution, delta_f])

    freqs_subset = freqs[(freqs >= (frequency - delta_f)) & (freqs <= (frequency + delta_f))]

    #corresponding fft bins
    lower_bin = np.nonzero(freqs >= freqs_subset[0])[0][0]
    upper_bin = np.nonzero(freqs >= freqs_subset[-1])[0][0] + 1 #+1 due to range stopping short of end, but need to include it
    num_bins = upper_bin - lower_bin
    if num_bins == 0:
        num_bins = 1
        upper_bin = lower_bin + num_bins

    #load samples
    timestamps, samples = load_gnuradio_samples(gnuradio_filename, return_full_timestamps=False)

    if num_seconds is not None:
        seconds = (timestamps[1] - timestamps[0])/np.timedelta64(1, 's')
        samples_per_second = len(samples)/(1.0*seconds)
        num_samples = int(samples_per_second * num_seconds)
        samples = samples[0:num_samples]
        num_seconds = num_samples/samples_per_second
        timestamps = [timestamps[0], timestamps[0] + num_seconds*np.timedelta64(1, 's')]

    #calculate spectrum at the desired subset
    fft = spectrum(timestamps, samples, n_fft=n_fft, start_bin=lower_bin, num_bins=num_bins, return_timestamps=False, save_memory=True)

    #change the timestamps to start and end in the center of first and last FFT windows
    time_resolution = (timestamps[-1] - timestamps[0])/(len(samples)*1.0)
    start_timestamp = timestamps[0] + n_fft/2.0*time_resolution
    end_timestamp = timestamps[0] + (len(fft)*n_fft - n_fft/2.0)*time_resolution
    timestamps = np.array([start_timestamp, end_timestamp])

    return freqs[lower_bin:upper_bin], fft, timestamps

def spectrum(timestamps, iq_samples, start_bin=0, num_bins=None, n_fft=8192, return_timestamps=True, save_memory=False):
    """
    Calculate FFT spectrum of a given signal. Blackman-Harris windowing is applied.

    Parameters
    ----------
    timestamps :
        Timestamps, for changing the rate of the timestamps to the rate of the
        FFT windows
    iq_samples :
        Raw IQ samples as obtained from GNU Radio
    start_bin : optional
        Start bin for the bins we want to keep
    num_bins : optional
        Number of bins we want to keep. If None, will return all bins
    n_fft : optional
        Size of FFT
    return_timestamps : optional
        Whether to return resampled timestamps
    save_memory : optional
        If False, will run the FFT transform in one go on the full data array
        and then take the subset specified by start_bin and num_bins, and
        otherwise take the FFT window by window and reduce to the subset for
        each window. The same end result is achieved, but the former is faster
        and the latter is less memory-extensive.
    Returns
    -------
    timestamps : optional
        Resampled timestamps, if return_timestamps is set to True
    ret_spectrum :
        FFT spectrum at specified bins
    """
    if num_bins is None:
        num_bins = n_fft

    w = blackman(n_fft)

    num_frames = len(iq_samples)/n_fft
    end_bin = start_bin + num_bins

    if not save_memory:
        ret_spectrum = fftshift(fft(np.reshape(iq_samples[0:num_frames*n_fft], (num_frames, n_fft)), n=n_fft), axes=1)[:, start_bin:end_bin]
    else:
        ret_spectrum = np.zeros((num_frames, num_bins), dtype=np.complex64)
        for i in np.arange(0, num_frames):
            fft_res = np.fft.fftshift(np.fft.fft(w*iq_samples[i*n_fft:(i+1)*n_fft]))
            ret_spectrum[i,:] = fft_res[start_bin:end_bin]

    if return_timestamps:
        timestamps = timestamps[n_fft/2:-1:n_fft]
        timestamps = timestamps[0:ret_spectrum.shape[0]]
        return timestamps, ret_spectrum
    else:
        return ret_spectrum

def get_center_frequency_and_samplerate(gnuradio_filename):
    """
    Get center frequency and sample rate in Hz from GNU Radio header.

    Assumed to be detached header, with header file in a .hdr file with gnuradio_filename as the filename without .hdr.

    Parameters
    ----------
    gnuradio_filename : string
        Filename for the IQ samples. Headerfile assumed to be gnuradio_filename.hdr
    Returns
    -------
    f_c :
        Center frequency in Hz
    f_s :
        Sample rate in Hz
    """

    header = load_gnuradio_header(gnuradio_filename + '.hdr')
    f_c = pmt.to_double(header['rx_freq'])
    f_s = header['rx_rate']
    return f_c, f_s

def frequency_axis(f_c, f_s, n_fft):
    """
    Calculate the frequency axis for an FFT.

    Parameters
    ----------
    f_c :
        Center frequency
    f_s :
        Samplerate
    n_fft :
        Number of FFT bins
    Returns
    -------
    fft_axis :
        Frequency axis for the corresponding (shifted) discrete fourier transform
    """
    return np.fft.fftshift(np.fft.fftfreq(n_fft))*f_s + f_c
