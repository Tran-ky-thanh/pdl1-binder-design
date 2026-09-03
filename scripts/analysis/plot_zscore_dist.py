import os
PROJECT = os.environ.get("PDL1_PROJECT", os.path.expanduser("~/protein_designs/pdl1_bench"))
OUTDIR = os.environ.get("PDL1_OUT", "build")
os.makedirs(OUTDIR, exist_ok=True)
import csv, statistics as st, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
SC=PROJECT+"/fold_screen/complex_run/scoring"
OUTW=OUTDIR
rows=[]
for r in csv.DictReader(open(SC+"/final_ranking.csv")):
    if r["composite"]:
        rows.append((int(r["rank"]),r["cand"],float(r["composite"])))
rows.sort(key=lambda x:x[0])
comp=[c for _,_,c in rows]
n=len(rows); med=st.median(comp); mean=st.mean(comp)
fig,ax=plt.subplots(1,2,figsize=(15,6))
# histogram
ax[0].hist(comp,bins=14,color="#3b7a57",edgecolor="w",alpha=.9)
ax[0].axvline(med,ls="--",c="#c33",lw=1.5,label="median %.2f"%med)
ax[0].axvline(0,ls=":",c="#888",lw=1,label="z=0 (pool mean)")
ax[0].set_xlabel("final composite z-score  (4·z(ipSAE) + 1·z(sc_DockQ))",fontsize=11)
ax[0].set_ylabel("số designs",fontsize=11)
ax[0].set_title("Phân bố z-score cuối (n=%d, 3-model)"%n,fontweight="bold",fontsize=12)
ax[0].legend(fontsize=9); ax[0].grid(alpha=.25)
# ranked
xs=[r for r,_,_ in rows]
colors=["#e8792b" if r<=10 else ("#4878b0" if r<=30 else "#b9c2cc") for r,_,_ in rows]
ax[1].bar(xs,comp,color=colors,width=.85)
ax[1].axhline(0,c="#888",lw=.8)
ax[1].axvline(10.5,ls="--",c="#e8792b",lw=1,alpha=.7); ax[1].text(10.6,max(comp)*.9,"top 10",color="#e8792b",fontsize=8)
ax[1].axvline(30.5,ls="--",c="#4878b0",lw=1,alpha=.7); ax[1].text(30.6,max(comp)*.9,"top 30",color="#4878b0",fontsize=8)
for r,c,v in rows[:5]: ax[1].annotate(c.replace("cand",""),(r,v),fontsize=7,ha="center",va="bottom")
# mark cand02117
for r,c,v in rows:
    if c=="cand02117": ax[1].annotate("02117",(r,v),fontsize=7,color="#a00",ha="center",va="top")
ax[1].set_xlabel("rank",fontsize=11); ax[1].set_ylabel("composite z-score",fontsize=11)
ax[1].set_title("Xếp hạng (cam=top10, xanh=top30, xám=còn lại)",fontweight="bold",fontsize=12)
ax[1].grid(alpha=.25,axis="y")
fig.suptitle("PD-L1 binders: final z-score distribution — Boltz-2 + ESMFold2-Full + Protenix",fontsize=13,fontweight="bold")
fig.tight_layout(rect=[0,0,1,.96])
fig.savefig(SC+"/final_zscore_dist.png",dpi=145); fig.savefig(OUTW+"/final_zscore_dist.png",dpi=145)
print("n=%d med=%.2f mean=%.2f min=%.2f max=%.2f"%(n,med,mean,min(comp),max(comp)))
q=lambda p: sorted(comp)[min(n-1,int(p*n))]
print("Q1=%.2f Q3=%.2f"%(q(.25),q(.75)))
print("gaps (rank: z, delta to next):")
for i in range(len(rows)-1):
    d=rows[i][2]-rows[i+1][2]
    if d>0.9: print("  after rank %d (%s z=%.2f): gap %.2f"%(rows[i][0],rows[i][1],rows[i][2],d))
print("saved final_zscore_dist.png")
