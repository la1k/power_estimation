"""
Probabilistic approach for antenna pattern extraction from a morse signal.
Decides whether a sample is a signal or not based on the noise floor
probability distribution.  The technique is described in
https://www.la1k.no/?p=2962.
"""

from scipy.ndimage.morphology import binary_erosion
import numpy as np
from scipy.stats import normaltest, foldnorm

def probablistic_signal_template_extraction(beacon_signal, non_beacon_signal):
    """
    Characterize the noise floor, and use its probability distribution to decide
    which samples of the beacon signal most likely are non-pauses.

    Parameters
    ----------
    beacon_signal : complex
        Beacon signal (as obtained from an FFT bin)
    non_beacon_signal : complex
        Non-beacon signal (as obtained from an FFT bin sufficiently far away in
        frequency from the beacon signal, but not so far that the noise floor
        no longer has the same distribution)

    Returns
    -------
    high_values : boolean
        Signal template, with 1 corresponding to high values, 0 to low values
    """

    #run test for normality on segments of the noise floor
    window_length = 5000
    window_pos = np.arange(window_length, len(non_beacon_signal) - window_length)

    normaltest_pvalues = np.zeros(len(window_pos))
    for i, start_samp in enumerate(window_pos):
        subset = non_beacon_signal[start_samp-window_length/2:start_samp+window_length/2]
        _, normaltest_pvalues[i] = normaltest(np.real(subset))

    #select samples that within a stretch of normal distributed samples, with critical value of 0.3, use this to characterize the noise floor distribution
    window_pos_normal = window_pos[normaltest > 0.3]
    normal_samples = non_beacon_signal[window_pos_normal]
    mean_noisefloor = np.mean(np.real(normal_samples))
    std_noisefloor = np.std(np.real(normal_samples))

    #get folded normal distribution (distribution of signal magnitudes of the noise)
    c = 0
    loc = 0
    S = std_noisefloor
    distrib = foldnorm(c=0, scale=S, loc=loc)

    #calculate probability that the noise floor probability density produces a
    #more extreme value than each observed sample
    noise_floor_prob = 1-distrib.cdf(np.abs(beacon_signal))

    #signal threshold (implicit threshold: probabilities lower than the
    #granularity of the datatype). Might have to be adjusted for more low-level
    #signals.
    high_values = noise_floor_prob == 0.0
    high_values = binary_erosion(high_values)

    return high_values
