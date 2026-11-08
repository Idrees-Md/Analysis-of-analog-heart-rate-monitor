import json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, welch, spectrogram
D="data/"
FILES={"rest_01":"ppg_rest_01.csv","rest_02":"ppg_rest_02.csv","rest_03":"ppg_rest_03.csv",
"strapped_01":"ppg_strapped_01.csv","strapped_02":"ppg_strapped_02.csv",
"noise_01":"ppg_noise_01.csv","paper_01":"ppg_paper_01.csv"}
FINGER=["rest_01","rest_02","rest_03","strapped_01","strapped_02"]; NOFINGER=["noise_01","paper_01"]
fs=200; W=30*fs
rng=np.random.default_rng(20261108)
INK,SEC,GRID,GREY="#0b0b0b","#52514e","#e4e3df","#8a8985"
BLUE,ORANGE,AQUA,YELLOW="#2a78d6","#eb6834","#1baf7a","#eda100"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.edgecolor":SEC,"axes.labelcolor":INK,"xtick.color":SEC,"ytick.color":SEC,
 "axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":0.8,"axes.grid":True,"grid.color":GRID,"grid.linewidth":0.6,"axes.axisbelow":True,
 "figure.facecolor":"white","axes.facecolor":"white","savefig.dpi":200,"savefig.bbox":"tight","axes.titlesize":10,"axes.titleweight":"bold","axes.titlelocation":"left"})
def load(n):
    d=pd.read_csv(D+FILES[n]); ts=d.timestamp_us.values.astype(np.int64); v=d.value.values.astype(float)
    t=(ts-ts[0])/1e6; tu=np.arange(0,t[-1],1/fs); x=np.interp(tu,t,v)
    return x, np.diff(ts), v
DATA={}; LATE={}; RAWSTD={}; RAWMEAN={}
for n in FILES:
    x,dt,v=load(n); DATA[n]=x; LATE[n]=int((dt>7500).sum()); RAWSTD[n]=float(v.std()); RAWMEAN[n]=float(v.mean())
DET={n:x-np.polyval(np.polyfit(np.arange(len(x))/fs,x,1),np.arange(len(x))/fs) for n,x in DATA.items()}
b,a=butter(3,[0.6/100,5/100],'band'); b2,a2=butter(3,[0.6/100,4/100],'band')
R=json.load(open("results/sim_results.json"))
def fmt(ax,xl=None,yl=None):
    if xl: ax.set_xlabel(xl)
    if yl: ax.set_ylabel(yl)
def label(ax,x,y,s,c,**k): ax.text(x,y,s,color=c,fontsize=8.5,fontweight="bold",va="center",**k)
