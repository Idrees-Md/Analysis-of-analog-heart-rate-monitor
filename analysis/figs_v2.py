exec(open("analysis/common_setup.py").read())
EST=["FFT peak (0.7-3.5 Hz)","FFT peak (50-110 bpm)","Autocorrelation","Harmonic product spectrum"]
def bands(x):
    n=len(x); X=np.abs(np.fft.rfft((x-x.mean())*np.hanning(n)))**2; f=np.fft.rfftfreq(n,1/fs); tot=X.sum()
    return [X[(f>=lo)&(f<hi)].sum()/tot*100 for lo,hi in [(0,.5),(.5,4),(4,10),(10,40),(40,100.1)]]
order=["rest_01","rest_02","rest_03","strapped_01","strapped_02","noise_01","paper_01"]
BE={n:bands(DET[n]) for n in order}
def inband(seg): return float(np.std(filtfilt(b2,a2,seg-seg.mean())))
def synth(hr,dur,r):
    t=np.arange(0,dur,1/fs); sig=np.zeros_like(t); T=60/hr; tb=-r.uniform(0,T); beats=[]
    while tb<dur+T:
        Tb=T*(1+0.03*r.standard_normal())
        for amp,c,w in [(1.0,0.20,0.07),(0.4,0.50,0.10)]: sig+=amp*np.exp(-0.5*((t-(tb+c*Tb))/(w*Tb))**2)
        beats.append(tb+0.2*Tb); tb+=Tb
    return sig-sig.mean(), np.array(beats)
def savefig(fig,name): fig.savefig("results/figs/"+name); plt.close(fig)
# spectral peak of rest_03 for honest titles
x3=DET["rest_03"]; Xf=np.abs(np.fft.rfft(x3*np.hanning(len(x3)))); ff=np.fft.rfftfreq(len(x3),1/fs); m=(ff>=2)&(ff<=4); pk3=ff[m][np.argmax(Xf[m])]
# ---- Fig 1
fig,ax=plt.subplots(2,1,figsize=(7,4.4),sharex=True); t10=np.arange(0,10,1/fs)
for a_,n,c,ttl in [(ax[0],"noise_01",GREY,"No finger (sensor covered)"),(ax[1],"rest_03",BLUE,"Fingertip at rest")]:
    a_.plot(t10,DET[n][40*fs:50*fs],color=c,lw=0.7); a_.set_ylim(-300,300); a_.set_ylabel("ADC counts\n(detrended)")
    a_.text(0.01,0.93,ttl,transform=a_.transAxes,fontsize=9,fontweight="bold",color=INK,va="top")
ax[1].set_xlabel("Time (s)")
ratio=RAWSTD["rest_03"]/RAWSTD["noise_01"]
fig.suptitle(f"A finger makes the signal {ratio:.0f}x larger, but it oscillates near {pk3:.1f} Hz ({pk3*60:.0f} bpm)",x=0.125,y=0.985,ha="left",fontsize=10,fontweight="bold")
fig.text(0.125,-0.01,"Source: noise_01 and rest_03 (HW-827 + XIAO ESP32-C6, 200 Hz), 10 s each, same vertical scale. The hand count during rest_03 was about 78 bpm.",fontsize=7.5,color=SEC)
savefig(fig,"fig1_time_domain.png")
# ---- Fig 2
fig,ax=plt.subplots(figsize=(7,4.2)); ax.axvspan(0.6,4,color=BLUE,alpha=0.08,lw=0)
for n,c in [("noise_01",GREY),("rest_03",BLUE),("strapped_01",ORANGE)]:
    f,P=welch(DET[n],fs,window="hann",nperseg=2048,noverlap=1024); mk=f>=0.3; ax.loglog(f[mk],P[mk],color=c,lw=1.1)
ax.set_xlim(0.3,100); ax.set_ylim(1e-3,3e3)
ax.text(0.62,6e-3,"pulse band 0.6-4 Hz",fontsize=8,color=SEC)
label(ax,0.35,0.35,"no finger",GREY,ha="left"); label(ax,1.25,1700,"fingertip (rest_03)",BLUE,ha="left"); label(ax,1.25,520,"strapped (strapped_01)",ORANGE,ha="left")
ax.annotate("50 Hz mains",xy=(50,1.2),xytext=(50,0.12),fontsize=8,color=SEC,ha="center",arrowprops=dict(arrowstyle="-",color=SEC,lw=0.6))
fmt(ax,"Frequency (Hz)","Power spectral density (counts$^2$/Hz)")
ax.set_title("Finger-on power sits at 10-20 Hz, far above the 0.6-4 Hz pulse band",pad=8)
fig.text(0.125,-0.02,"Welch estimate, Hann window, 2048-sample segments, 50% overlap, linear detrend. Source: 120 s recordings at 200 Hz.",fontsize=7.5,color=SEC)
savefig(fig,"fig2_psd.png")
# ---- Fig 3
fig,ax=plt.subplots(figsize=(7,3.9))
cols=["#c9c8c3",BLUE,"#a9a8a3","#8a8985","#d9d8d3"]; names=["<0.5 Hz","0.5-4 Hz (pulse band)","4-10 Hz","10-40 Hz","40-100 Hz"]
for i,n in enumerate(order):
    y=len(order)-1-i; left=0
    for k,v in enumerate(BE[n]):
        ax.barh(y,v,left=left,color=cols[k],edgecolor="white",linewidth=1.5,height=0.62)
        if k==1: ax.text(left+v+1.2,y,f"{v:.1f}%",va="center",fontsize=8,color=BLUE,fontweight="bold",bbox=dict(fc="white",ec="none",pad=1.2))
        left+=v
ax.set_yticks(range(len(order))); ax.set_yticklabels([f"{n}  (std {RAWSTD[n]:.0f})" for n in order][::-1],fontsize=8.5)
ax.set_xlim(0,100); ax.set_xlabel("Share of signal energy (%)"); ax.grid(axis="y",visible=False)
ax.axhline(1.5,color=SEC,lw=0.6,ls=":"); ax.text(100,1.5,"below: no finger / no reflection",ha="right",va="center",fontsize=7,color=SEC,bbox=dict(fc="white",ec="none",pad=0.5))
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=cols[k],label=names[k]) for k in range(5)],ncol=5,loc="upper center",bbox_to_anchor=(0.45,1.14),frameon=False,fontsize=7.5,handlelength=1,columnspacing=1)
ax.set_title(f"Only {min(BE[n][1] for n in order[:5]):.1f}-{max(BE[n][1] for n in order[:5]):.1f}% of finger-on energy lies in the pulse band",pad=30)
fig.text(0.125,-0.02,"Hann-windowed FFT power, detrended. Std in ADC counts. No-finger files hold very little total energy (std 7-9 counts).",fontsize=7.5,color=SEC)
savefig(fig,"fig3_band_energy.png")
# ---- Fig 4
fig,axs=plt.subplots(1,2,figsize=(7.4,3.7),sharey=True); S={}
for n in ["rest_03","strapped_01"]:
    f,t,Sx=spectrogram(DET[n],fs,window="hann",nperseg=2000,noverlap=1800); S[n]=(f,t,10*np.log10(Sx+1e-12))
allv=np.concatenate([S[n][2][S[n][0]<6].ravel() for n in S]); vmin,vmax=np.percentile(allv,5),np.percentile(allv,99.5)
for a_,n,ttl in zip(axs,S,["rest_03 (fingertip, rest)","strapped_01 (strapped fingertip)"]):
    f,t,Z=S[n]; mk=f<=6; im=a_.pcolormesh(t,f[mk],Z[mk],cmap="Blues",vmin=vmin,vmax=vmax,shading="auto")
    for y in (50/60,110/60): a_.axhline(y,color=ORANGE,lw=0.9,ls="--")
    a_.set_title(ttl,fontsize=9,pad=6); a_.set_xlabel("Time (s)"); a_.grid(False)
axs[0].set_ylabel("Frequency (Hz)"); axs[0].text(8,1.33,"resting heart rate\n50-110 bpm",color=ORANGE,fontsize=7.5,va="center",fontweight="bold",bbox=dict(fc="white",ec="none",alpha=0.7,pad=1))
sec=axs[1].secondary_yaxis("right",functions=(lambda f:f*60,lambda b:b/60)); sec.set_ylabel("bpm")
cb=fig.colorbar(im,ax=axs,pad=0.12,shrink=0.85); cb.set_label("Power (dB)")
fig.suptitle("The strongest steady ridge sits near 3 Hz (about 180-210 bpm), outside the resting-rate band",x=0.06,y=1.0,ha="left",fontsize=10,fontweight="bold")
savefig(fig,"fig4_spectrogram.png")
# ---- Fig 5
R=json.load(open("results/sim_results.json")); snr=np.array(R["snrs"]); Nn=R["trials"]
def wilson(p,n,z=1.96):
    d=1+z*z/n; c=(p+z*z/(2*n))/d; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/d; return c-h,c+h
CE=dict(zip(EST,[BLUE,ORANGE,AQUA,YELLOW])); MK=dict(zip(EST,["o","s","^","D"]))
def cross90(y):
    y=np.array(y)
    for i in range(1,len(y)):
        if y[i-1]<0.9<=y[i]: return float(snr[i-1]+(snr[i]-snr[i-1])*(0.9-y[i-1])/(y[i]-y[i-1]))
    return float("nan")
X90={k:cross90(R["success"][k]) for k in EST}
fig,(axA,axB)=plt.subplots(2,1,figsize=(7,6.2),sharex=True,gridspec_kw={"height_ratios":[1.15,1]})
for k in EST:
    p=np.array(R["success"][k])*100; lo,hi=zip(*[wilson(v/100,Nn) for v in p]); fc="white" if k.startswith("FFT peak (0.7") else CE[k]
    axA.plot(snr,p,color=CE[k],lw=1.5,marker=MK[k],ms=5,mfc=fc,mew=1.3,label=k); axA.vlines(snr,np.array(lo)*100,np.array(hi)*100,color=CE[k],lw=0.8,alpha=0.6)
    axB.plot(snr,R["mae"][k],color=CE[k],lw=1.5,marker=MK[k],ms=5,mfc=fc,mew=1.3)
axA.axhline(90,color=SEC,lw=0.8,ls="--"); axA.text(-24.6,92.5,"90% success",fontsize=8,color=SEC)
axA.set_ylim(-3,103); axA.set_ylabel("Correct within 3 bpm (%)"); axA.legend(loc="lower right",frameon=False,fontsize=8)
axB.set_yscale("log"); axB.set_ylabel("Mean absolute error (bpm)"); axB.set_xlabel("Pulse-to-noise ratio in the 0.6-4 Hz band (dB)"); axB.set_xticks(snr)
axA.set_title(f"A pulse must reach {X90['FFT peak (50-110 bpm)']:.0f} to {X90['Autocorrelation']:+.0f} dB of the in-band noise for 90% success",pad=8)
fig.text(0.125,0.005,f"{Nn} trials per point (vertical bars: 95% Wilson interval). Synthetic two-Gaussian pulse, 55-105 bpm, 30 s windows, added to real finger-on noise windows.",fontsize=7.5,color=SEC)
savefig(fig,"fig5_simulation_accuracy.png")
# ---- Fig 6
H=R["hist"]; bins=np.array(H["bins"]); ctr=(bins[:-1]+bins[1:])/2
fig,ax=plt.subplots(figsize=(7,3.7)); w=1.9
ax.bar(ctr-w/2,np.array(H["finger"])*100,width=w,color=BLUE,label=f"finger on ({H['nw']} windows)"); ax.bar(ctr+w/2,np.array(H["nofinger"])*100,width=w,color=GREY,label=f"no finger ({H['nn']} windows)")
share=lambda h:sum(v for v,c in zip(h,ctr) if 60<=c<80)*100
ax.set_xlabel("Estimated heart rate (bpm); harmonic product spectrum, 50-110 bpm search"); ax.set_ylabel("Share of 30 s windows (%)"); ax.legend(frameon=False,fontsize=8,loc="upper right"); ax.set_xticks(bins[::2])
ax.set_title(f"Fed noise alone, the estimator answers 60-80 bpm in {share(H['finger']):.0f}% (finger) and {share(H['nofinger']):.0f}% (no finger) of windows",pad=8,fontsize=9.5)
fig.text(0.125,-0.03,f"Real recordings only, no synthetic pulse. Median {H['med_f']:.0f} bpm (finger) and {H['med_n']:.0f} bpm (no finger). 5 bpm bins.",fontsize=7.5,color=SEC)
savefig(fig,"fig6_noise_only_histogram.png")
# ---- Fig 7 (band-passed view, plus raw std for contrast)
base=DET["rest_03"][:W]; p,beats=synth(72,30,np.random.default_rng(7)); prms=inband(p); nrms=inband(base)
fig,axs=plt.subplots(3,1,figsize=(7,5.4),sharex=True)
for a_,s in zip(axs,[6,-6,-18]):
    x=base+p*(nrms*10**(s/20)/prms); y=filtfilt(b2,a2,x-x.mean())[2*fs:12*fs]
    a_.plot(np.arange(2*fs,12*fs)/fs,y,color=BLUE,lw=0.9); a_.set_ylim(-48,56); a_.set_ylabel("counts\n(0.6-4 Hz)")
    bb=beats[(beats>=2)&(beats<12)]; a_.plot(bb,np.full_like(bb,49),marker="v",ls="none",color=ORANGE,ms=5)
    a_.text(0.01,0.12,f"pulse-to-noise {s:+d} dB",transform=a_.transAxes,fontsize=9,fontweight="bold",va="bottom",color=INK,bbox=dict(fc="white",ec="none",alpha=0.8,pad=1))
axs[0].text(0.99,0.12,"orange marks: true beats (72 bpm)",transform=axs[0].transAxes,color=ORANGE,fontsize=8,fontweight="bold",ha="right",va="bottom",bbox=dict(fc="white",ec="none",alpha=0.8,pad=1))
axs[-1].set_xlabel("Time (s)")
fig.suptitle("Inside the pulse band, beats can be followed at +6 dB, are marginal at -6 dB and are lost at -18 dB",x=0.125,y=0.965,ha="left",fontsize=10,fontweight="bold")
fig.text(0.125,0.005,f"Noise: first 30 s of rest_03, band-passed 0.6-4 Hz (noise RMS {nrms:.1f} counts in band; raw std about {np.std(base):.0f}). In the raw trace the pulse is not visible at any of these levels.",fontsize=7.2,color=SEC)
savefig(fig,"fig7_synthetic_examples.png")
# ---- Fig 8
f=np.linspace(0.01,500,5000); Hf=lambda x:abs(np.sin(5*np.pi*x/1000)/(5*np.sin(np.pi*x/1000))); Hb=Hf(f)
fig,ax=plt.subplots(figsize=(7,3.7)); ax.plot(f,Hb,color=BLUE,lw=1.4)
for x in (4,17,50,100): ax.plot([x],[Hf(x)],marker="o",color=ORANGE,ms=4)
for x in (200,400): ax.plot([x],[0],marker="o",color=ORANGE,ms=4); ax.text(x,0.06,f"null\n{x} Hz",fontsize=8,ha="center",color=INK)
ax.text(255,0.80,"Gain at marked points:\n4 Hz     "+f"{Hf(4):.3f}\n17 Hz   {Hf(17):.3f}\n50 Hz   {Hf(50):.3f}\n100 Hz {Hf(100):.3f}",fontsize=8.5,color=INK,va="top",family="DejaVu Sans Mono")
ax.set_xlim(0,500); ax.set_ylim(0,1.12); fmt(ax,"Frequency at the 1 kHz sample rate (Hz)","Gain |H(f)|")
ax.set_title("A 5-sample average leaves the pulse band intact and nulls 200 and 400 Hz",pad=8)
fig.text(0.125,-0.03,"H(f) = sin(5 pi f/1000) / (5 sin(pi f/1000)). It attenuates other frequencies only mildly, so it is a weak anti-alias filter.",fontsize=7.5,color=SEC)
savefig(fig,"fig8_boxcar_response.png")
# ---- prominence + equivalent SNR + Fig 9 + tables
def prom(seg):
    y=filtfilt(b,a,seg-seg.mean()); n=len(y); Xs=np.abs(np.fft.rfft(y*np.hanning(n),n=4*n)); ff=np.fft.rfftfreq(4*n,1/fs)
    m=(ff>=50/60)&(ff<=110/60); mm=(ff>=0.5)&(ff<=4); return float(Xs[m].max()/np.median(Xs[mm]))
def wins(names): return [DET[n][s:s+W] for n in names for s in range(0,len(DET[n])-W,10*fs)]
fw,nw=wins(FINGER),wins(NOFINGER); pn_f=np.array([prom(s) for s in fw]); pn_n=np.array([prom(s) for s in nw])
SN=[-24,-18,-12,-6,0,6,12]; med=[];q1=[];q3=[]
for s_ in SN:
    v=[]
    for _ in range(100):
        seg=fw[rng.integers(len(fw))]; hr=rng.uniform(55,105); pp,_=synth(hr,30,rng); pp=pp[:W]; v.append(prom(seg+pp*(inband(seg)*10**(s_/20)/inband(pp))))
    med.append(np.median(v)); q1.append(np.percentile(v,25)); q3.append(np.percentile(v,75))
med=np.array(med); REAL={n:float(np.median([prom(DET[n][s:s+W]) for s in range(0,len(DET[n])-W,10*fs)])) for n in order}
def equiv(p):
    if p<=np.percentile(pn_f,95): return None
    if p>=med[-1]: return float("inf")
    return float(np.interp(p,med,SN))
fig,ax=plt.subplots(figsize=(7,4.2)); ax.fill_between(SN,q1,q3,color=BLUE,alpha=0.15,lw=0); ax.plot(SN,med,color=BLUE,lw=1.5,marker="o",ms=4)
ax.set_yscale("log"); ax.set_ylim(1.5,45)
p50,p95=np.median(pn_f),np.percentile(pn_f,95)
ax.axhline(p50,color=GREY,lw=1,ls="--"); ax.text(-3,p50*0.88,f"noise only, finger on: median {p50:.1f}",fontsize=8,color=SEC,va="top")
ax.axhline(p95,color=GREY,lw=1,ls=":"); ax.text(-3,p95*1.07,f"noise only, finger on: 95th percentile {p95:.1f}",fontsize=8,color=SEC,va="bottom")
ax.axhline(max(REAL.values()),color=ORANGE,lw=1.1,ls="--"); ax.text(-6.5,max(REAL.values())*0.93,f"highest real recording ({max(REAL,key=REAL.get)}): median {max(REAL.values()):.1f}",fontsize=8,color=ORANGE,fontweight="bold",va="top")
ax.set_xlabel("Simulated pulse-to-noise ratio (dB)"); ax.set_ylabel("Peak prominence (peak / median spectrum)"); ax.set_xticks(SN)
ax.set_title("Every real recording's peak prominence stays below the noise-only 95th percentile",pad=8,fontsize=9.5)
fig.text(0.125,-0.03,"Prominence = tallest FFT peak in 50-110 bpm over the median spectrum in 0.5-4 Hz, 30 s windows. Band = interquartile range, 100 trials per level. Log scale.",fontsize=7.5,color=SEC)
savefig(fig,"fig9_peak_prominence.png")
T=R["table"]; rows=[]
cond={"rest_01":"fingertip, rest (200 Hz direct)","rest_02":"fingertip, rest (1 kHz avg to 200 Hz)","rest_03":"fingertip, rest (hand count ~78 bpm)","strapped_01":"strapped fingertip","strapped_02":"strapped fingertip","noise_01":"no finger, covered","paper_01":"against paper, no reflection"}
for n in order:
    rows.append(dict(recording=n,condition=cond[n],samples=T[n]["rows"],late_samples=T[n]["late"],mean_counts=round(T[n]["mean"],1),std_counts=round(T[n]["std"],1),pct_below_0p5Hz=round(BE[n][0],1),pct_0p5_4Hz=round(BE[n][1],1),pct_4_10Hz=round(BE[n][2],1),pct_10_40Hz=round(BE[n][3],1),pct_40_100Hz=round(BE[n][4],1),top_peak_above_4Hz=round(T[n]["top_peak_above4"],1),inband_rms_counts=round(inband(DET[n]),2),median_prominence=round(REAL[n],2)))
pd.DataFrame(rows).to_csv("results/tables/table_recordings.csv",index=False)
pd.DataFrame({"snr_db":R["snrs"],**{k:R["success"][k] for k in EST}}).to_csv("results/tables/sim_success_rate.csv",index=False)
pd.DataFrame({"snr_db":R["snrs"],**{k:R["mae"][k] for k in EST}}).to_csv("results/tables/sim_mae_bpm.csv",index=False)
pd.DataFrame({"bin_start_bpm":H["bins"][:-1],"finger_on_share":H["finger"],"no_finger_share":H["nofinger"]}).to_csv("results/tables/noise_only_histogram.csv",index=False)
eq=[dict(recording=n,median_prominence=round(REAL[n],2),equivalent_snr_db=("not above noise-only 95th percentile" if equiv(REAL[n]) is None else round(equiv(REAL[n]),1))) for n in order]
pd.DataFrame(eq).to_csv("results/tables/equivalent_snr.csv",index=False)
nf_rms=float(np.median([inband(s) for s in fw])); no_rms=float(np.median([inband(s) for s in nw])); raw_std_f=float(np.median([np.std(s) for s in fw]))
json.dump(dict(X90=X90,nf_rms=nf_rms,no_rms=no_rms,raw_std_finger_median=raw_std_f,REAL=REAL,prom_noise_finger=[float(p50),float(p95)],prom_noise_nofinger=[float(np.median(pn_n)),float(np.percentile(pn_n,95))],med=[float(x) for x in med],eq=eq,Hf={str(x):float(Hf(x)) for x in (4,17,50,100)},pk3=float(pk3)),open("results/tables/calc_numbers.json","w"),indent=1)
print("X90",{k:round(v,2) for k,v in X90.items()}); print("inband rms finger/no finger",round(nf_rms,2),round(no_rms,2),"raw std finger median",round(raw_std_f,1)); print("prom noise finger p50/p95",round(p50,2),round(p95,2),"| no finger",round(np.median(pn_n),2),round(np.percentile(pn_n,95),2)); print(pd.DataFrame(eq).to_string()); print("pk3",pk3)
