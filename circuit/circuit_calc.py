import numpy as np, json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
f=np.logspace(-2,2.5,3000); s=2j*np.pi*f
def par(a,b): return a*b/(a+b)
R1,C1,Rf,Cf=33e3,10e-6,330e3,100e-9
def stage(Rin_extra=0):
    Zin=R1+Rin_extra+1/(s*C1); Zf=par(Rf,1/(s*Cf)); return -Zf/Zin
out={"corners":{"HP_Hz":round(1/(2*np.pi*R1*C1),2),"LP_Hz":round(1/(2*np.pi*Rf*Cf),2),"stage_gain":Rf/R1}}
fig,ax=plt.subplots(figsize=(7,3.8),dpi=200)
cols={1.0:"#2a78d6",0.5:"#1baf7a",0.1:"#eb6834"}
for k in (1.0,0.5,0.1):
    Rs=k*(1-k)*10e3
    H=stage()*k*stage(Rs); db=20*np.log10(abs(H)); pk=db.max(); band=f[db>=pk-3]
    ax.semilogx(f,db,color=cols[k],lw=2,label=f"gain preset wiper at {int(k*100)}%")
    out[f"k={k}"]={"peak_dB":round(float(pk),1),"peak_x":round(float(10**(pk/20))),
      "f_lo_3dB":round(float(band[0]),2),"f_hi_3dB":round(float(band[-1]),2),
      "dB_at_0.1Hz":round(float(np.interp(0.1,f,db)),1),"dB_at_15Hz":round(float(np.interp(15,f,db)),1),
      "dB_at_50Hz":round(float(np.interp(50,f,db)),1),"dB_at_1.5Hz":round(float(np.interp(1.5,f,db)),1),
      "rel_15Hz_vs_peak":round(float(np.interp(15,f,db)-pk),1),"rel_50Hz_vs_peak":round(float(np.interp(50,f,db)-pk),1)}
ax.axvspan(0.6,4,color="#8a8985",alpha=.15,label="pulse band 0.6–4 Hz")
ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Gain (dB)"); ax.set_ylim(-60,50); ax.grid(alpha=.25,which="both")
ax.set_title("Designed analog front end: calculated response (two inverting band-pass stages)",fontsize=9,loc="left")
ax.legend(frameon=False,fontsize=8,loc="lower center"); ax.spines[["top","right"]].set_visible(False)
fig.tight_layout(); fig.savefig("results/circuit/fig_afe_response.png")
out["hysteresis_mV"]=round(3.5*10/(10+680)*1000)
out["divider"]={"ratio":round(100/394,3),"Rthev_k":round(294*100/394,1)}
out["led_mA"]=round((5-2.1-0.2)/470*1000,1)
out["threshold_range_V"]=[round(5*6.8/23.6,2),round(5*16.8/23.6,2)]
out["Vref"]=2.5
json.dump(out,open("results/circuit/circuit_numbers.json","w"),indent=1); print(json.dumps(out,indent=1))
