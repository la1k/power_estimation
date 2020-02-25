Scripts, Python modules and GNU Radio companion examples for measuring the
reception power of an antenna during rotation, and for doing the necessary
off-line analysis to produce antenna diagrams.

![Example pattern](example-pattern.png)

The code here is described in a blog post, https://www.la1k.no/2018/05/16/using-gnu-radio-and-hamlib-to-calculate-antenna-diagrams/.

Use `measure.py` in measurement subdirectory to sample IQ data while
simultaneously logging antenna angles to file. See `README.md` in that folder
for more details. Let the script run continuously while turning the rotor.

Use `combine_samples_and_angles.py` to combine samples/power estimates produced
using GNU Radio with the above angle measurements.  See docstring of
`combine_samples_and_angles.py` for more details.

For time-varying, periodic morse signals, see docstring of `cw_extraction.py` for extracting of the high signal levels only.
Its behavior is described in:

* https://www.la1k.no/2018/06/27/extracting-the-antenna-pattern-from-a-beacon-signal-pt-i-initial-investigations-using-autocorrelation/
* https://www.la1k.no/2018/07/04/extracting-the-antenna-pattern-from-a-beacon-signal-pt-ii-measurements-of-la2shf/
* https://www.la1k.no/2018/08/01/extracting-the-antenna-pattern-from-a-beacon-signal-pt-iii-frequency-analysis-and-code/

See also `probabilistic_signal_extraction.py` for a basic probabilistic approach. Its behavior is described in https://www.la1k.no/2018/08/08/extracting-the-antenna-pattern-from-a-beacon-signal-pt-iv-probabilistic-approach/.
