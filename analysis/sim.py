import numpy as np, pandas as pd, json
from scipy.signal import butter, filtfilt, welch
D="data/"
FILES={"rest_01":"ppg_rest_01.csv","rest_02":"ppg_rest_02.csv","rest_03":"ppg_rest_03.csv",
"strapped_01":"ppg_strapped_01.csv","strapped_02":"ppg_strapped_02.csv",
"noise_01":"ppg_noise_01.csv","paper_01":"ppg_paper_01.csv"}
FINGER=["rest_01","rest_02","rest_03","strapped_01","strapped_02"]; NOFINGER=["noise_01","paper_01"]
fs=200
rng=np.random.default_rng(20261107)
def load(name):
    d=pd.read_csv(D+FILES[name]); ts=d.timestamp_us.values.astype(np.int64); v=d.value.values.astype(float)
    t=(ts-ts[0])/1e6; dt=np.diff(ts)
    tu=np.arange(0,t[-1],1/fs); x=np.interp(tu,t,v)
    raw=dict(rows=len(v),dur=float(t[-1]),late=int((dt>7500).sum()),mean=float(v.mean()),std=float(v.std()),vmin=float(v.min()),vmax=float(v.max()))
    x=x-np.polyval(np.polyfit(tu,x,1),tu)
    return x,raw
data={n:load(n) for n in FILES}
b,a=butter(3,[0.6/(fs/2),5.0/(fs/2)],'band')
b2,a2=butter(3,[0.6/(fs/2),4.0/(fs/2)],'band')
def inband_rms(seg): return float(np.std(filtfilt(b2,a2,seg-seg.mean())))
# ---- estimators on a 30 s raw segment
def spec(y):
    n=len(y); X=np.abs(np.fft.rfft(y*np.hanning(n),n=4*n)); f=np.fft.rfftfreq(4*n,1/fs); return f,X
def est_fft_naive(y):
    f,X=spec(y); m=(f>=0.7)&(f<=3.5); return 60*f[m][np.argmax(X[m])]
def est_fft_range(y):
    f,X=spec(y); m=(f>=50/60)&(f<=110/60); return 60*f[m][np.argmax(X[m])]
def est_hps(y):
    f,X=spec(y); n3=len(X)//3; H=X[:n3]*X[::2][:n3]*X[::3][:n3]; ff=f[:n3]
    m=(ff>=50/60)&(ff<=110/60); return 60*ff[m][np.argmax(H[m])]
def est_acf(y):
    n=len(y); y=y-y.mean(); F=np.fft.rfft(y,2*n); ac=np.fft.irfft(np.abs(F)**2)[:n]; ac/=ac[0]
    lag=np.arange(n)/fs; m=(lag>=60/130)&(lag<=60/45)
    idx=np.where(m)[0]; seg=ac[idx]
    pk=[i for i in range(1,len(seg)-1) if seg[i]>seg[i-1] and seg[i]>=seg[i+1]]
    if not pk: return np.nan
    j=max(pk,key=lambda i:seg[i]); return 60/lag[idx[j]]
EST={"FFT peak (0.7-3.5 Hz)":est_fft_naive,"FFT peak (50-110 bpm)":est_fft_range,"Autocorrelation":est_acf,"Harmonic product spectrum":est_hps}
def run_est(seg):
    y=filtfilt(b,a,seg-seg.mean()); return {k:f(y) for k,f in EST.items()}
# ---- synthetic PPG (two-Gaussian beat model, 3% beat-to-beat interval variability)
def synth_ppg(hr,dur,rng):
    t=np.arange(0,dur,1/fs); sig=np.zeros_like(t); T=60/hr; tb=-rng.uniform(0,T)
    while tb<dur+T:
        Tb=T*(1+0.03*rng.standard_normal())
        for amp,c,w in [(1.0,0.20,0.07),(0.4,0.50,0.10)]:
            sig+=amp*np.exp(-0.5*((t-(tb+c*Tb))/(w*Tb))**2)
        tb+=Tb
    return sig-sig.mean()
# ---- noise bank: real 30 s finger-on windows (step 10 s)
W=30*fs
def windows(names):
    out=[]
    for n in names:
        x=data[n][0]
        for s in range(0,len(x)-W,10*fs): out.append((n,x[s:s+W]))
    return out
fw=windows(FINGER); nw=windows(NOFINGER)
print("finger windows",len(fw),"no-finger windows",len(nw))
med_noise_rms=float(np.median([inband_rms(s) for _,s in fw])); med_floor_rms=float(np.median([inband_rms(s) for _,s in nw]))
print("median in-band (0.6-4 Hz) noise RMS counts: finger",round(med_noise_rms,2),"no finger",round(med_floor_rms,2))
SNRS=[-24,-18,-12,-6,0,6,12]; N=200; TOL=3.0
success={k:[] for k in EST}; mae={k:[] for k in EST}
for snr in SNRS:
    ok={k:0 for k in EST}; err={k:[] for k in EST}
    for _ in range(N):
        _,seg=fw[rng.integers(len(fw))]
        hr=rng.uniform(55,105); p=synth_ppg(hr,30,rng)[:W]
        prms=float(np.std(filtfilt(b2,a2,p))); nrms=inband_rms(seg)
        x=seg+p*(nrms*10**(snr/20)/prms)
        r=run_est(x)
        for k,v in r.items():
            e=abs(v-hr); err[k].append(e); ok[k]+= e<=TOL
    for k in EST: success[k].append(ok[k]/N); mae[k].append(float(np.nanmean(err[k])))
    print(snr,{k[:12]:round(ok[k]/N,2) for k in EST})
# ---- pure-noise estimates (real windows, no synthetic pulse)
BINS=list(range(50,115,5))
def hist(wins,fn):
    v=np.array([fn(filtfilt(b,a,s-s.mean())) for _,s in wins]); h,_=np.histogram(v,bins=BINS); return (h/len(v)).tolist(), float(np.median(v)), float(np.std(v))
h_f,m_f,s_f=hist(fw,est_hps); h_n,m_n,s_n=hist(nw,est_hps)
print("HPS on finger-on windows: median",round(m_f,1),"sd",round(s_f,1),"| no-finger: median",round(m_n,1),"sd",round(s_n,1))
# ---- measured characterization table
def bands(x):
    n=len(x); X=np.abs(np.fft.rfft((x-x.mean())*np.hanning(n)))**2; f=np.fft.rfftfreq(n,1/fs); tot=X.sum()
    return [float(X[(f>=lo)&(f<hi)].sum()/tot*100) for lo,hi in [(0,.5),(.5,4),(4,10),(10,40),(40,100.1)]]
table={}
for n in FILES:
    x,raw=data[n]; f_,X_=spec(x)
    m=f_>4; top=float(f_[m][np.argmax(X_[m])])
    table[n]=dict(raw,bands=bands(x),top_peak_above4=top)
# ---- PSD for chart (log-spaced bins)
def psd_rows(name):
    x=data[name][0]; f,P=welch(x,fs,window='hann',nperseg=2048,noverlap=1024)
    edges=np.logspace(np.log10(0.5),np.log10(99),41); out=[]
    for i in range(40):
        m=(f>=edges[i])&(f<edges[i+1])
        if m.any(): out.append([float(np.sqrt(edges[i]*edges[i+1])),float(np.log10(P[m].mean()))])
    return out
psd={n:psd_rows(n) for n in ["noise_01","rest_03","strapped_01"]}
res=dict(snrs=SNRS,success=success,mae=mae,trials=N,tol=TOL,med_noise_rms=med_noise_rms,med_floor_rms=med_floor_rms,
 hist=dict(bins=BINS,finger=h_f,nofinger=h_n,med_f=m_f,sd_f=s_f,med_n=m_n,sd_n=s_n,nw=len(fw),nn=len(nw)),table=table,psd=psd)
json.dump(res,open("results/sim_results.json","w"),indent=1)
