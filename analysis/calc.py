import json, numpy as np, pandas as pd
from scipy.signal import butter, filtfilt
exec(open("analysis/common_setup.py").read())
C=json.load(open("results/tables/calc_numbers.json")); R=json.load(open("results/sim_results.json")); T=pd.read_csv("results/tables/table_recordings.csv")
Hf=lambda x:abs(np.sin(5*np.pi*x/1000)/(5*np.sin(np.pi*x/1000)))
ff=np.linspace(0.1,500,200000); f3=float(ff[np.argmax(Hf(ff)<=1/np.sqrt(2))])
def spec(y):
    n=len(y); X=np.abs(np.fft.rfft(y*np.hanning(n),n=4*n)); f=np.fft.rfftfreq(4*n,1/fs); return f,X
def hps(y):
    f,X=spec(y); n3=len(X)//3; H=X[:n3]*X[::2][:n3]*X[::3][:n3]; g=f[:n3]; m=(g>=50/60)&(g<=110/60); return 60*g[m][np.argmax(H[m])]
hr={}
for n in ["rest_01","rest_02","rest_03","strapped_01","strapped_02","noise_01","paper_01"]:
    v=[hps(filtfilt(b,a,DET[n][s:s+W]-DET[n][s:s+W].mean())) for s in range(0,len(DET[n])-W,10*fs)]; hr[n]=(float(np.median(v)),float(np.std(v)),len(v))
nf,no,rawstd=C["nf_rms"],C["no_rms"],C["raw_std_finger_median"]
X90=C["X90"]; dc=float(np.mean([RAWMEAN[n] for n in FINGER])); lsb=3.3/4095*1000
L=[]
w=L.append
w("# Calculations for the paper\n")
w("All numbers are computed from the seven recordings (HW-827 + XIAO ESP32-C6, 200 Hz, 120 s each) and from the simulation in `sim.py`. Scripts are in the zip. Nothing here uses MAX30102 data: none was collected.\n")
w("## 1. Sampling and frequency resolution\n")
w(f"- Sampling rate fs = 200 Hz, so the Nyquist limit is fs/2 = 100 Hz.\n- Recording length T = 120 s gives N = 24,000 samples. Bin spacing = fs/N = 1/T = {1/120:.5f} Hz = {60/120:.2f} bpm.\n- A 30 s analysis window gives {1/30:.4f} Hz = {60/30:.1f} bpm per bin. To resolve 1 bpm (0.0167 Hz) you need T >= 60 s.\n- Zero-padding by 4 (used in the estimators) only interpolates between bins. It does not improve resolution.\n- Mains hum at 50 Hz is below Nyquist (100 Hz), so it appears at 50 Hz and can be filtered. Sampling at 100 Hz would put it exactly at the Nyquist limit.\n")
w("## 2. ADC scale (nominal)\n")
w(f"- 12-bit ADC: 4096 levels. If full scale is taken as 3.3 V, one count = 3.3/4095 = {lsb:.3f} mV. The real range depends on the ESP32-C6 attenuation setting, so treat millivolt values as nominal.\n- Typical finger-on mean level: {dc:.0f} counts (about {dc*lsb/1000:.2f} V nominal).\n- No-finger noise: std {RAWSTD['noise_01']:.1f} counts (noise_01) and {RAWSTD['paper_01']:.1f} counts (paper_01), about {RAWSTD['noise_01']*lsb:.1f} mV nominal. This is the electronics floor.\n")
w("## 3. Timing quality\n")
w("| Recording | Samples | Samples later than 7.5 ms | Share |\n|---|---|---|---|")
for _,r in T.iterrows(): w(f"| {r.recording} | {r.samples:,} | {r.late_samples} | {r.late_samples/r.samples*100:.2f}% |")
w("\nOnly rest_01 (first firmware, 200 Hz direct sampling) has no late samples. The later firmware (1 kHz sampling averaged to 200 Hz) shows 59 to 60 late samples (one 8 ms interval about every 2 s). They were repaired by interpolation onto an even 200 Hz grid before analysis.\n")
w("## 4. The 5-sample average (1 kHz to 200 Hz)\n")
w(f"- Gain of a 5-point moving average at fs = 1 kHz: |H(f)| = sin(5 pi f/1000) / (5 sin(pi f/1000)).\n- Random (white) noise falls by sqrt(5) = {np.sqrt(5):.2f}x ({20*np.log10(np.sqrt(5)):.1f} dB) in amplitude.\n- Gain at 4 Hz = {Hf(4):.4f}, 17 Hz = {Hf(17):.4f}, 50 Hz = {Hf(50):.3f}, 100 Hz = {Hf(100):.3f}. The -3 dB point is at about {f3:.0f} Hz.\n- Nulls at 200 Hz and 400 Hz. Other frequencies are only mildly attenuated: gain at 150 Hz = {Hf(150):.3f}, at 250 Hz = {Hf(250):.3f}.\n- After averaging, the 200 Hz output can still carry aliases: a component at 183 Hz (gain {Hf(183):.3f}) folds to 17 Hz, and one at 217 Hz (gain {Hf(217):.3f}) also folds to 17 Hz. The average is a weak anti-alias filter, so it is not a substitute for analog filtering ahead of the ADC.\n")
w("## 5. Measured noise characterization\n")
w("| Recording | Condition | Std (counts) | In-band RMS 0.6-4 Hz (counts) | Energy 0.5-4 Hz | Energy 10-40 Hz | Top peak above 4 Hz (Hz) | Median peak prominence |\n|---|---|---|---|---|---|---|---|")
for _,r in T.iterrows(): w(f"| {r.recording} | {r.condition} | {r.std_counts} | {r.inband_rms_counts} | {r.pct_0p5_4Hz}% | {r.pct_10_40Hz}% | {r.top_peak_above_4Hz} | {r.median_prominence} |")
w(f"\n- Median in-band (0.6-4 Hz) noise RMS over 30 s windows: finger on = {nf:.2f} counts, no finger = {no:.2f} counts. Ratio {nf/no:.2f} ({20*np.log10(nf/no):.1f} dB).\n- Median raw std of finger-on windows (detrended): {rawstd:.1f} counts. The in-band noise is {nf/rawstd*100:.1f}% of that in amplitude ({(nf/rawstd)**2*100:.1f}% in power), consistent with the 1.0 to 2.7% band energy above.\n- The strongest narrow peaks above 4 Hz sit at 12 to 18 Hz and move between recordings (12.0, 14.3, 16.2, 16.6, 17.6 Hz). Their source is not identified.\n")
w("## 6. Simulation: how strong must the pulse be?\n")
w("Method: a two-Gaussian pulse (systolic and diastolic components, 3% beat-to-beat interval variability, true rate drawn uniformly from 55 to 105 bpm) is scaled so its 0.6-4 Hz RMS equals the stated ratio times the in-band RMS of a real finger-on noise window, then added to that window. 200 trials per level, 30 s windows. Success means the estimate is within 3 bpm of the truth.\n")
w("| Pulse-to-noise (dB) | "+" | ".join(R["success"].keys())+" |\n|---|"+"---|"*len(R["success"]))
for i,s in enumerate(R["snrs"]): w(f"| {s:+d} | "+" | ".join(f"{R['success'][k][i]*100:.0f}%" for k in R["success"])+" |")
w("\nPulse-to-noise ratio at which each estimator reaches 90% success (linear interpolation between simulated levels), and the in-band pulse amplitude that implies for your finger-on noise:\n")
w("| Estimator | 90% point (dB) | Pulse RMS needed (counts, 0.6-4 Hz) | As % of the typical level | Nominal mV | Raw-trace SNR (dB) |\n|---|---|---|---|---|---|")
for k,v in X90.items():
    amp=nf*10**(v/20); w(f"| {k} | {v:+.1f} | {amp:.1f} | {amp/dc*100:.2f}% | {amp*lsb:.1f} | {20*np.log10(amp/rawstd):.1f} |")
w(f"\nThe raw-trace SNR column divides the needed in-band pulse RMS by the raw std ({rawstd:.1f} counts). It shows why the pulse is invisible in the time-domain trace even when an estimator can recover it.\n")
w("## 7. What an estimator reports on noise alone\n")
H=R["hist"]; ctr=[(a+b)/2 for a,b in zip(H["bins"][:-1],H["bins"][1:])]
sh=lambda h:sum(v for v,c in zip(h,ctr) if 60<=c<80)*100
w(f"- Harmonic product spectrum on real windows with no synthetic pulse: finger on ({H['nw']} windows) median {H['med_f']:.0f} bpm, sd {H['sd_f']:.1f}; no finger ({H['nn']} windows) median {H['med_n']:.0f} bpm, sd {H['sd_n']:.1f}.\n- Share of windows reporting 60 to 80 bpm: {sh(H['finger']):.0f}% (finger) and {sh(H['nofinger']):.0f}% (no finger).\n- So a reported heart rate of 60 to 75 bpm from these recordings carries no evidence of a pulse. The earlier 'about 69 bpm' reading is what the estimator returns on noise.\n")
w("## 8. Real recordings against the noise-only level\n")
w(f"- Peak prominence = tallest FFT peak in 50-110 bpm divided by the median spectrum in 0.5-4 Hz (30 s windows).\n- Noise only, finger on: median {C['prom_noise_finger'][0]:.2f}, 95th percentile {C['prom_noise_finger'][1]:.2f}. No finger: median {C['prom_noise_nofinger'][0]:.2f}, 95th percentile {C['prom_noise_nofinger'][1]:.2f}.\n- Median prominence of the real recordings: "+", ".join(f"{n} {v:.2f}" for n,v in C['REAL'].items())+f". All are below the finger-on 95th percentile ({C['prom_noise_finger'][1]:.2f}), so none shows a peak distinguishable from noise.\n- Simulated median prominence rises from {C['med'][0]:.1f} at -24 dB to {C['med'][3]:.1f} at -6 dB and {C['med'][4]:.1f} at 0 dB.\n")
w("## 9. The one reference check\n")
w(f"- A hand count during rest_03 gave 39 beats in 30 s = 78 bpm (self-reported; plus or minus one beat is plus or minus 2 bpm).\n- Harmonic product spectrum on rest_03, 30 s windows: median {hr['rest_03'][0]:.0f} bpm, sd {hr['rest_03'][1]:.1f} over {hr['rest_03'][2]} windows. It does not reproduce 78 bpm.\n- The same estimator on the other recordings (median bpm, sd): "+"; ".join(f"{n} {m:.0f} ({s:.1f})" for n,(m,s,_) in hr.items())+". The no-finger recordings give the same range as the finger-on ones.\n")
w("## 10. Definitions\n")
w("- Pulse-to-noise ratio (dB) = 20 log10(pulse RMS / noise RMS), both measured in the 0.6-4 Hz band after a 3rd-order Butterworth band-pass applied forwards and backwards.\n- Success rate interval: Wilson 95% interval for a binomial proportion with n = 200.\n- Peak prominence: see section 8.\n- Estimators: (1) strongest FFT peak 0.7-3.5 Hz; (2) strongest FFT peak 50-110 bpm; (3) autocorrelation peak for 45-130 bpm; (4) harmonic product spectrum with 3 harmonics, 50-110 bpm. All use a Hann window, 4x zero-padding (1, 2, 4) and a 0.6-5 Hz band-pass.\n")
w("## Limits to state in the paper\n")
w("- The finger-on recordings may contain a weak real pulse. The simulation adds a synthetic pulse on top of them, so the in-band noise is slightly overestimated.\n- Only one recording (rest_03) has a reference, and it is a hand count.\n- The synthetic pulse is a generic two-Gaussian model, not fitted to this sensor.\n- Results apply to this sensor, this placement and these recordings. They say nothing about the MAX30102.\n")
open("results/calculations.md","w").write("\n".join(L))
print("\n".join(L)[:6000])
