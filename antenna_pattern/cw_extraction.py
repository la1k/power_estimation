"""
Functions for extracting the signal level from a morse signal, for antenna pattern measurements.

See extract_signal_level_from_periodic_cw_signal() as an entry point method.
"""

from frequency_analysis import interpolate_at_frequency
import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage import binary_erosion

def extract_signal_level_from_periodic_cw_signal(filename, frequency, n_fft=8192):
    """
    Extract the high signal level from a periodic CW signal, skipping over the pauses.

    Parameters
    ----------
    filename : string
        GNU Radio metafile filename (assuming detached header)
    frequency :
        Frequency of the CW signal in Hz
    n_fft : optional
        Length of FFT

    Returns
    -------
    high_stamps :
        Timestamps for samples corresponding to high level of the signal
    high_samples :
        Samples corresponding to high level
    """

    #get signal exactly at the given frequency
    stamps, samples = interpolate_at_frequency(filename, frequency, n_fft=n_fft)

    time_resolution = (stamps[-1] - stamps[0])/(1.0*len(samples))
    full_stamps = np.arange(0, len(samples))*time_resolution + stamps[0]

    #divide into repeating sequences
    sequences = divide_periodic_cw_signal_in_segments(stamps, samples)

    #rectify sequence offsets
    rectification_lags, sequences = rectify_sequences_by_crosscorrelation(sequences)

    #get signal template
    template = construct_signal_template(sequences, rectification_lags)[0:len(samples)]
    template = binary_erosion(template > 0)

    #apply signal template
    high_samples = samples[template > 0]
    high_stamps = full_stamps[template > 0]

    return high_stamps, high_samples

def divide_periodic_cw_signal_in_segments(stamps, samples):
    """
    Assuming that the input signal is periodic, find the period and divide the
    signal into a two-dimensional array with the sequences along the first
    axis.

    Parameters
    ----------
    stamps : tuple
        Specifies the first and last timestamp of the signal, as obtained from
        load_gnuradio_samples with return_full_timestamps set to False (and
        after potential transformation by the FFT).
    samples : complex, array_like
        Complex signal to be decomposed into sequences.

    Returns
    -------
    segment_image : ndarray
        The signal divided into its periodic sequences as an
        num_periods x samples_in_period matrix.
    """

    seconds = (stamps - stamps[0])/np.timedelta64(1, 's')

    #assume that the CW signal is machine-generated and extremely periodic:
    #find period using autocorrelation. Need to use logarithms in order to make
    #the entire pattern comparable against itself
    ma_len = 100
    corr, lags = apply_autocorrelation(np.log(np.abs(samples)) - moving_average(np.log(np.abs(samples)), ma_len))
    corr = corr[lags >= 0]
    lags = lags[lags >= 0]

    #assume first autocorrelation peak to appear at lag of more than one
    #second, and that the peak is very clear and will appear cleanly when we
    #sort the values
    sorted_correlations = np.flip(np.argsort(corr), axis=-1)
    time_resolution = seconds[1]/(1.0*len(samples))
    peak_thresh = 1/time_resolution
    signal_period = int(sorted_correlations[sorted_correlations > peak_thresh][0])

    #divide the signal into periodic signal segments
    num_segments = int(np.ceil(len(samples)/(1.0*signal_period)))
    segment_image = np.zeros((num_segments, signal_period))
    for i in np.arange(0, num_segments):
        sequence = samples[i*signal_period:(i+1)*signal_period]
        segment_image[i,0:len(sequence)] = sequence
    return segment_image

def rectify_sequences_by_crosscorrelation(segment_image):
    """
    Use crosscorrelation to rectify offsets between the sequences obtained from
    divide_periodic_cw_signal_in_segments() so that the sequences match.

    Parameters
    ----------
    segment_image : ndarray
        Sequence image as obtained from divide_periodic_cw_signal_in_segments().

    Returns
    -------
    corr_lags : ndarray
        The timing offsets of each sequence as compared to the first sequence.
    segment_image_rectified : ndarray
        Rectified sequence image.
    """

    def normalize(seq):
        """
        Normalization for application of cross-correlation.
        """
        seq = np.nan_to_num(np.log(seq))
        return seq - np.mean(np.invert(np.isnan(seq)))

    #apply cross-correlation for all sequences against the first sequence
    #in order to find the lags
    segment_image_rectified = np.ones(segment_image.shape)
    corr_lags = np.zeros(len(segment_image_rectified))
    _, lags = apply_autocorrelation(segment_image[0, :])
    for i in np.arange(1, segment_image.shape[0]-1):
        seq_1_orig = segment_image[0,:]
        seq_2_orig = segment_image[i,:]

        seq_1 = normalize(seq_1_orig)
        seq_2 = normalize(seq_2_orig)

        correlation = np.correlate(seq_1, seq_2, 'full')
        maxind = np.argmax(correlation)
        corr_lags[i] = lags[maxind]

    #use the lags to rectify the offset of each sequence
    for i, lag in enumerate(corr_lags):
        seg_length = len(segment_image[0,:])
        old_x = np.arange(0, seg_length)
        rectified_x = old_x + int(lag)

        valid_x = (rectified_x >= 0) & (rectified_x < seg_length)
        segment_image_rectified[i, rectified_x[valid_x]] = segment_image[i, valid_x]
    return corr_lags, segment_image_rectified

def construct_signal_template(segment_image, rectification_lags=None):
    """
    Estimate where the signal is high and low.

    Parameters
    ----------
    segment_image : ndarray
        Sequence image, as obtained from divide_periodic_cw_signal_in_segments() (or rectified in rectify_sequences_by_crosscorrelation()).
    rectification_lags : ndarray, optional
        Offsets of the sequences in segment_image with respect to the original signal, as obtained from rectify_sequences_by_crosscorrelation().

    Returns
    -------
    samples_template : boolean, ndarray
        Template, where 1 corresponds to signal high and 0 corresponds to signal low in the underlying morse signal of the original signal.
    """

    if rectification_lags is None:
        rectification_lags = np.zeros(segment_image.shape[0])

    num_segments = segment_image.shape[0]
    signal_period = segment_image.shape[1]
    sequence_template = np.max(segment_image, axis=0) > np.mean(segment_image)
    samples_template = np.zeros(num_segments*len(sequence_template))

    for i in np.arange(0, num_segments):
        start = i*signal_period
        end = (i+1)*signal_period

        #estimate approximate signal level throughout this specific sequence
        sequence = segment_image[i, :]
        high = sequence_template > 0
        x = np.arange(0, len(sequence))
        signal_level = moving_average(interp1d(x[high], sequence[high])(x[high]), 50)

        #normalize the sequence by the signal level approximation (but don't normalize what we know already is background)
        rect_signal = sequence.copy()
        rect_signal[sequence_template > 0] = sequence[sequence_template > 0]/signal_level

        #estimate signal and background means
        signal_mean = np.mean(rect_signal[sequence_template])
        background_mean = np.mean(rect_signal[np.invert(sequence_template)])
        threshold = (signal_mean + background_mean)/2

        #improve on the general template by mean-thresholding of this sequence
        rect_template = sequence_template.copy()
        local_threshold = rect_signal > threshold
        rect_template = local_threshold & sequence_template

        #place local template at the correct position as according to the
        #cross-correlated offset
        offset = int(rectification_lags[i])
        samples_template[start-offset:end-offset] = rect_template

    return samples_template

def moving_average(x, window_len):
    """
    Calculate central moving average of a signal.

    From scipy-cookbook, modified to return the central moving average.

    Parameters
    ----------
    x :
        Signal of which to calculate moving average
    window_len :
        Length of moving average window
    Returns
    -------
    y :
        Central moving average of the input signal
    """
    w = np.ones(window_len, 'd')
    s = np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    y = np.convolve(w/w.sum(), s, mode='valid')

    #get offset from being the CMA
    offset = 0
    if ((window_len % 2) == 0):
        offset = int(window_len/2)
    else:
        offset = int((window_len - 1)/2)

    return y[offset:-offset+1]

def apply_autocorrelation(signal):
    """
    Calculate autocorrelation of a signal.

    Parameters
    ----------
    signal : ndarray
        Signal to be autocorrelated.

    Returns
    -------
    autocorr : ndarray
        Autocorrelation
    lags : ndarray
        Corresponding lags for each autocorrelation value
    """

    maxlag = len(signal)-1
    lags = np.array(range(-maxlag, maxlag+1))

    autocorr = np.correlate(signal, signal, 'full')
    return autocorr, lags
