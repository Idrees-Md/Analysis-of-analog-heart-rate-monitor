# Calculations for the paper

All numbers are computed from the seven recordings (HW-827 + XIAO ESP32-C6, 200 Hz, 120 s each) and from the simulation in `sim.py`. Scripts are in the zip. Nothing here uses MAX30102 data: none was collected.

## 1. Sampling and frequency resolution

- Sampling rate fs = 200 Hz, so the Nyquist limit is fs/2 = 100 Hz.
- Recording length T = 120 s gives N = 24,000 samples. Bin spacing = fs/N = 1/T = 0.00833 Hz = 0.50 bpm.
- A 30 s analysis window gives 0.0333 Hz = 2.0 bpm per bin. To resolve 1 bpm (0.0167 Hz) you need T >= 60 s.
- Zero-padding by 4 (used in the estimators) only interpolates between bins. It does not improve resolution.
- Mains hum at 50 Hz is below Nyquist (100 Hz), so it appears at 50 Hz and can be filtered. Sampling at 100 Hz would put it exactly at the Nyquist limit.

## 2. ADC scale (nominal)

- 12-bit ADC: 4096 levels. If full scale is taken as 3.3 V, one count = 3.3/4095 = 0.806 mV. The real range depends on the ESP32-C6 attenuation setting, so treat millivolt values as nominal.
- Typical finger-on mean level: 1951 counts (about 1.57 V nominal).
- No-finger noise: std 6.8 counts (noise_01) and 8.6 counts (paper_01), about 5.5 mV nominal. This is the electronics floor.

## 3. Timing quality

| Recording | Samples | Samples later than 7.5 ms | Share |
|---|---|---|---|
| rest_01 | 24,001 | 0 | 0.00% |
| rest_02 | 23,973 | 60 | 0.25% |
| rest_03 | 23,966 | 59 | 0.25% |
| strapped_01 | 23,966 | 59 | 0.25% |
| strapped_02 | 23,966 | 59 | 0.25% |
| noise_01 | 23,973 | 60 | 0.25% |
| paper_01 | 23,965 | 60 | 0.25% |

Only rest_01 (first firmware, 200 Hz direct sampling) has no late samples. The later firmware (1 kHz sampling averaged to 200 Hz) shows 59 to 60 late samples (one 8 ms interval about every 2 s). They were repaired by interpolation onto an even 200 Hz grid before analysis.

## 4. The 5-sample average (1 kHz to 200 Hz)

- Gain of a 5-point moving average at fs = 1 kHz: |H(f)| = sin(5 pi f/1000) / (5 sin(pi f/1000)).
- Random (white) noise falls by sqrt(5) = 2.24x (7.0 dB) in amplitude.
- Gain at 4 Hz = 0.9994, 17 Hz = 0.9886, 50 Hz = 0.904, 100 Hz = 0.647. The -3 dB point is at about 90 Hz.
- Nulls at 200 Hz and 400 Hz. Other frequencies are only mildly attenuated: gain at 150 Hz = 0.312, at 250 Hz = 0.200.
- After averaging, the 200 Hz output can still carry aliases: a component at 183 Hz (gain 0.097) folds to 17 Hz, and one at 217 Hz (gain 0.084) also folds to 17 Hz. The average is a weak anti-alias filter, so it is not a substitute for analog filtering ahead of the ADC.

## 5. Measured noise characterization

| Recording | Condition | Std (counts) | In-band RMS 0.6-4 Hz (counts) | Energy 0.5-4 Hz | Energy 10-40 Hz | Top peak above 4 Hz (Hz) | Median peak prominence |
|---|---|---|---|---|---|---|---|
| rest_01 | fingertip, rest (200 Hz direct) | 82.7 | 10.99 | 2.7% | 50.8% | 16.2 | 3.29 |
| rest_02 | fingertip, rest (1 kHz avg to 200 Hz) | 96.9 | 11.16 | 1.4% | 74.5% | 12.0 | 2.13 |
| rest_03 | fingertip, rest (hand count ~78 bpm) | 80.3 | 6.97 | 1.1% | 66.3% | 14.3 | 2.23 |
| strapped_01 | strapped fingertip | 87.7 | 7.83 | 1.0% | 84.2% | 16.6 | 2.15 |
| strapped_02 | strapped fingertip | 74.9 | 9.2 | 2.1% | 64.3% | 17.6 | 2.81 |
| noise_01 | no finger, covered | 6.8 | 3.28 | 36.0% | 9.3% | 4.8 | 2.08 |
| paper_01 | against paper, no reflection | 8.6 | 4.69 | 48.4% | 12.0% | 5.8 | 2.06 |

- Median in-band (0.6-4 Hz) noise RMS over 30 s windows: finger on = 8.66 counts, no finger = 4.02 counts. Ratio 2.15 (6.7 dB).
- Median raw std of finger-on windows (detrended): 64.5 counts. The in-band noise is 13.4% of that in amplitude (1.8% in power), consistent with the 1.0 to 2.7% band energy above.
- The strongest narrow peaks above 4 Hz sit at 12 to 18 Hz and move between recordings (12.0, 14.3, 16.2, 16.6, 17.6 Hz). Their source is not identified.

## 6. Simulation: how strong must the pulse be?

Method: a two-Gaussian pulse (systolic and diastolic components, 3% beat-to-beat interval variability, true rate drawn uniformly from 55 to 105 bpm) is scaled so its 0.6-4 Hz RMS equals the stated ratio times the in-band RMS of a real finger-on noise window, then added to that window. 200 trials per level, 30 s windows. Success means the estimate is within 3 bpm of the truth.

| Pulse-to-noise (dB) | FFT peak (0.7-3.5 Hz) | FFT peak (50-110 bpm) | Autocorrelation | Harmonic product spectrum |
|---|---|---|---|---|
| -24 | 2% | 19% | 10% | 18% |
| -18 | 4% | 40% | 12% | 29% |
| -12 | 8% | 72% | 14% | 56% |
| -6 | 52% | 95% | 27% | 90% |
| +0 | 99% | 100% | 60% | 100% |
| +6 | 100% | 100% | 95% | 100% |
| +12 | 100% | 100% | 100% | 100% |

Pulse-to-noise ratio at which each estimator reaches 90% success (linear interpolation between simulated levels), and the in-band pulse amplitude that implies for your finger-on noise:

| Estimator | 90% point (dB) | Pulse RMS needed (counts, 0.6-4 Hz) | As % of the typical level | Nominal mV | Raw-trace SNR (dB) |
|---|---|---|---|---|---|
| FFT peak (0.7-3.5 Hz) | -1.1 | 7.6 | 0.39% | 6.1 | -18.6 |
| FFT peak (50-110 bpm) | -7.3 | 3.7 | 0.19% | 3.0 | -24.7 |
| Autocorrelation | +5.2 | 15.7 | 0.80% | 12.6 | -12.3 |
| Harmonic product spectrum | -6.1 | 4.3 | 0.22% | 3.5 | -23.5 |

The raw-trace SNR column divides the needed in-band pulse RMS by the raw std (64.5 counts). It shows why the pulse is invisible in the time-domain trace even when an estimator can recover it.

## 7. What an estimator reports on noise alone

- Harmonic product spectrum on real windows with no synthetic pulse: finger on (49 windows) median 69 bpm, sd 7.5; no finger (19 windows) median 68 bpm, sd 11.8.
- Share of windows reporting 60 to 80 bpm: 82% (finger) and 74% (no finger).
- So a reported heart rate of 60 to 75 bpm from these recordings carries no evidence of a pulse. The earlier 'about 69 bpm' reading is what the estimator returns on noise.

## 8. Real recordings against the noise-only level

- Peak prominence = tallest FFT peak in 50-110 bpm divided by the median spectrum in 0.5-4 Hz (30 s windows).
- Noise only, finger on: median 2.38, 95th percentile 4.64. No finger: median 2.07, 95th percentile 2.45.
- Median prominence of the real recordings: rest_01 3.29, rest_02 2.13, rest_03 2.23, strapped_01 2.15, strapped_02 2.81, noise_01 2.08, paper_01 2.06. All are below the finger-on 95th percentile (4.64), so none shows a peak distinguishable from noise.
- Simulated median prominence rises from 2.5 at -24 dB to 5.7 at -6 dB and 12.0 at 0 dB.

## 9. The one reference check

- A hand count during rest_03 gave 39 beats in 30 s = 78 bpm (self-reported; plus or minus one beat is plus or minus 2 bpm).
- Harmonic product spectrum on rest_03, 30 s windows: median 63 bpm, sd 12.0 over 10 windows. It does not reproduce 78 bpm.
- The same estimator on the other recordings (median bpm, sd): rest_01 68 (2.5); rest_02 70 (9.3); rest_03 63 (12.0); strapped_01 68 (3.7); strapped_02 73 (3.0); noise_01 68 (4.9); paper_01 62 (16.2). The no-finger recordings give the same range as the finger-on ones.

## 10. Definitions

- Pulse-to-noise ratio (dB) = 20 log10(pulse RMS / noise RMS), both measured in the 0.6-4 Hz band after a 3rd-order Butterworth band-pass applied forwards and backwards.
- Success rate interval: Wilson 95% interval for a binomial proportion with n = 200.
- Peak prominence: see section 8.
- Estimators: (1) strongest FFT peak 0.7-3.5 Hz; (2) strongest FFT peak 50-110 bpm; (3) autocorrelation peak for 45-130 bpm; (4) harmonic product spectrum with 3 harmonics, 50-110 bpm. All use a Hann window, 4x zero-padding (1, 2, 4) and a 0.6-5 Hz band-pass.

## Limits to state in the paper

- The finger-on recordings may contain a weak real pulse. The simulation adds a synthetic pulse on top of them, so the in-band noise is slightly overestimated.
- Only one recording (rest_03) has a reference, and it is a hand count.
- The synthetic pulse is a generic two-Gaussian model, not fitted to this sensor.
- Results apply to this sensor, this placement and these recordings. They say nothing about the MAX30102.
