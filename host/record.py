"""Record the serial stream from the XIAO ESP32-C6 to a CSV.
Edit PORT / FILE / SECS. Close the Arduino Serial Monitor first.
Find the port:  python -m serial.tools.list_ports
"""
import serial, time
PORT = "COM5"; FILE = "ppg_rest_02.csv"; SECS = 120
ser = serial.Serial(PORT, 921600, timeout=1)
with open(FILE, "w") as f:
    f.write("timestamp_us,value\n"); end = time.time() + SECS
    while time.time() < end:
        line = ser.readline().decode(errors="ignore").strip(); parts = line.split(",")
        if len(parts) != 2: continue
        try: int(parts[0]); float(parts[1])
        except ValueError: continue
        f.write(line + "\n")
print("saved", FILE)
