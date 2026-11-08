"""Heart-rate estimate from a PPG CSV (timestamp_us,value).
Usage: python hr.py ppg_rest_03.csv [min_bpm max_bpm]
Prints a bpm estimate per 15 s window plus a confidence score (0-1).
Windows with confidence < 0.30 are marked UNRELIABLE - do not trust those."""
import sys, numpy as np, pandas as pd
from scipy.signal import butter, filtfilt, find_peaks

FS = 200.0            # output rate after resampling
WIN, STEP = 15, 5     # seconds
MIN_CONF = 0.30       # heuristic threshold on autocorrelation height

path = sys.argv[1]
lo_bpm = float(sys.argv[2]) if len(sys.argv) > 2 else 45
hi_bpm = float(sys.argv[3]) if len(sys.argv) > 3 else 130

d = pd.read_csv(path)
t = (d.iloc[:, 0].values - d.iloc[0, 0]) / 1e6
v = d.iloc[:, 1].values.astype(float)
tu = np.arange(0, t[-1], 1 / FS)
x = np.interp(tu, t, v)                       # even 200 Hz grid (fixes dropped ticks)
x -= np.polyval(np.polyfit(tu, x, 1), tu)     # remove drift

# pulse band only: throws away the 10-40 Hz junk before anything else
b, a = butter(3, [0.6 / (FS / 2), 4.0 / (FS / 2)], btype="band")
y = filtfilt(b, a, x)

def estimate(seg):
    seg = seg - seg.mean()
    ac = np.correlate(seg, seg, "full")[len(seg) - 1:]
    ac /= ac[0]
    lag = np.arange(len(seg)) / FS
    m = (lag >= 60 / hi_bpm) & (lag <= 60 / lo_bpm)
    p, _ = find_peaks(ac[m])
    if len(p) == 0:
        return np.nan, 0.0
    best = p[np.argmax(ac[m][p])]
    return 60 / lag[m][best], float(ac[m][best])

rows = []
for s in np.arange(0, tu[-1] - WIN, STEP):
    seg = y[int(s * FS): int((s + WIN) * FS)]
    bpm, conf = estimate(seg)
    rows.append((s, bpm, conf))

print(f"{path}: search range {lo_bpm:.0f}-{hi_bpm:.0f} bpm")
print(" window(s)    bpm   confidence")
good = []
for s, bpm, conf in rows:
    flag = "" if conf >= MIN_CONF else "  UNRELIABLE"
    if conf >= MIN_CONF: good.append(bpm)
    print(f" {s:4.0f}-{s+WIN:<4.0f}  {bpm:6.1f}   {conf:.2f}{flag}")
if good:
    print(f"\nMedian of {len(good)}/{len(rows)} reliable windows: {np.median(good):.1f} bpm"
          f" (spread {np.std(good):.1f})")
else:
    print(f"\nNo reliable windows (0/{len(rows)}). Do not report a heart rate from this file.")
