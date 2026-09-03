import csv, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
R="/home/thanh/protein_designs/pdl1_bench"
rows=list(csv.DictReader(open(R+"/colabfold/cf_eval.csv")))
def short(n):
    n=n.replace("mythos_preview_single_target_pdl1_","M.st_").replace("mythos_preview_multi_target_pdl1_","M.mt_")
    n=n.replace("opus_4_8_multi_target_pdl1_","O.mt_")
    return n
x=np.array([float(r["ds_consensus"]) for r in rows])
y1=np.array([float(r["cf_ipsae_min"]) for r in rows])
y2=np.array([float(r["cf_iptm"]) for r in rows])
binder=[r["wetlab_binder"]=="True" for r in rows]
lab=[short(r["full_name"]) for r in rows]
r1=np.corrcoef(x,y1)[0,1]; r2=np.corrcoef(x,y2)[0,1]
print("r_ipSAE=%.3f  r_ipTM=%.3f"%(r1,r2))

fig,ax=plt.subplots(1,2,figsize=(12.4,5.2))
GRN="#1a9850"; RED="#c1272d"
def panel(a,y,ylabel,title,r):
    a.plot([0,1],[0,1],"--",color="#444",lw=1.2,zorder=1,label="y = x (exact agreement)")
    for xi,yi,b,l in zip(x,y,binder,lab):
        if b: a.scatter(xi,yi,s=95,facecolor=GRN,edgecolor=GRN,zorder=3)
        else: a.scatter(xi,yi,s=95,facecolor="none",edgecolor=RED,linewidths=1.8,zorder=3)
        a.annotate(l,(xi,yi),fontsize=7,xytext=(4,3),textcoords="offset points",color="#333")
    a.set_xlim(0,1); a.set_ylim(0,1); a.set_aspect("equal")
    a.set_xlabel("Dataset consensus ipSAE_min  (10 cofolding models)",fontsize=10)
    a.set_ylabel(ylabel,fontsize=10)
    a.set_title("%s\nPearson r = %.3f"%(title,r),fontsize=11)
    a.grid(alpha=.25)
# legend proxies on left panel
from matplotlib.lines import Line2D
leg=[Line2D([0],[0],marker='o',color='w',markerfacecolor=GRN,markersize=9,label='wetlab binder'),
     Line2D([0],[0],marker='o',color='w',markerfacecolor='none',markeredgecolor=RED,markeredgewidth=1.8,markersize=9,label='wetlab non-binder'),
     Line2D([0],[0],ls='--',color='#444',label='y = x (exact agreement)')]
panel(ax[0],y1,"localColabFold ipSAE_min","localColabFold ipSAE_min",r1)
panel(ax[1],y2,"localColabFold ipTM","localColabFold ipTM",r2)
ax[0].legend(handles=leg,fontsize=8.5,loc="upper left",framealpha=.95)
fig.suptitle("Dataset consensus vs. independent localColabFold — 7 PD-L1 designs",fontsize=13)
fig.tight_layout(rect=[0,0,1,.95])
out="/mnt/c/Users/thanh/AppData/Local/Temp/claude/C--Users-thanh-Documents/9a8e1c59-098d-4338-a1cc-e43c6d31654e/scratchpad/colabfold_consensus.png"
fig.savefig(out,dpi=150,bbox_inches="tight")
fig.savefig(R+"/colabfold/consensus_vs_localcolabfold.png",dpi=150,bbox_inches="tight")
import os; print("saved",out,round(os.path.getsize(out)/1024),"KB")
