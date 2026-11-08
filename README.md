# Heart-rate sensing: analog front end vs. embedded DSP

Signal-quality study and simulation for a low-cost pulse (PPG) heart-rate monitor.
Two approaches are compared: an analog electronics chain (op-amp filtering, comparator, LED alerts)
and an embedded/IoT chain (ESP8266 / ESP32-C6, inter-beat-interval BPM, Telegram alerts).

## Honest scope (read first)

- **No validated heart-rate accuracy is claimed.** The HW-827 recordings in `data/` did not yield a
  trustworthy pulse (best attempt: 63 ± 12 bpm vs. a hand count of ~78 bpm).
- The recordings are used to **characterise noise**, and a **Monte-Carlo simulation** shows how strong the
  pulse must be before each estimator works, and what estimators report on noise alone.
- The analog circuit in `circuit/` is a **designed reference circuit** with *calculated* (not measured)
  response. It is not a trace of the breadboard prototype.
- No MAX30102 data is included. A digital sensor is future work.

## Layout

| Folder | Contents |
|---|---|
| `firmware/ppg_logger/` | XIAO ESP32-C6 sketch: 1 kHz timer sampling, 5-sample average, 200 Hz CSV over serial |
| `host/record.py` | Records the serial stream to CSV |
| `data/` | Seven real HW-827 recordings (`timestamp_us,value`, 200 Hz): rest x3, strapped x2, no-finger noise, paper |
| `analysis/` | `figs_v2.py` (figures + tables), `calc.py` (calculations), `sim.py` (simulation), `hr.py` (estimator with confidence score) |
| `circuit/` | Designed analog front end: calculation and schematic scripts |
| `results/` | Figures 1-9, tables, `calculations.md`, circuit schematic and response |

## Reproduce

```bash
pip install -r requirements.txt
python analysis/sim.py          # simulation  -> results/sim_results.json
python analysis/figs_v2.py      # figures 1-9 and tables
python analysis/calc.py         # results/calculations.md
python circuit/circuit_calc.py  # front-end response + numbers
python circuit/schem.py         # schematic
python analysis/hr.py data/ppg_rest_03.csv
```
Run all commands from the repository root.

## Hardware

- Seeed XIAO ESP32-C6, HW-827 analog pulse sensor: S -> D0 (GPIO0), + -> 3V3, - -> GND.
- Sampling: 1 kHz hardware timer, averaged to 200 Hz (acts as a weak anti-alias filter; see Fig. 8).
- Designed front end (5 V): two inverting band-pass stages (0.48-4.8 Hz each), LM358, 200k/47k/47k + 100k divider to ESP8266 A0.

## Key figures

| | |
|---|---|
| Fig 1-4 | Time domain, PSD, band energy, spectrogram of the real recordings |
| Fig 5-7, 9 | Simulation: estimator success vs. pulse strength, noise-only histogram, synthetic examples, peak prominence |
| Fig 8 | Frequency response of the 5-sample moving average |
| `results/circuit/` | Designed circuit schematic and calculated response |

## License

MIT (see `LICENSE`).
