"""Plot the full 31-slot budget and patient-resampling sensitivity."""
from pathlib import Path
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import pandas as pd
import numpy as np

root = Path(sys.argv[1])
data = pd.read_csv(root / "confirmation.tsv", sep="\t")
diff = pd.read_csv(root / "differences.tsv", sep="\t")
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"svg.fonttype":"none","svg.hashsalt":"panel-31"})
fig = plt.figure(figsize=(8.0,6.1), layout="constrained")
grid = fig.add_gridspec(2,2,height_ratios=[2.2,1.05])
ink="#173f4a"; pale="#b5cbd0"; orange="#a24420"
labels=["Joint only","RNA alone","Protein alone","Joint overall","RNA alone","Protein alone"]
ys=[5.8,4.8,3.8,2,1,0]
for col,(ep,title) in enumerate([("RNA","A   Erasmus tumor RNA"),("protein","B   Huashan tumor protein")]):
    ax=fig.add_subplot(grid[0,col])
    ordered=[]
    for design,joint in [("primary","joint_only"),("exploratory","joint_overall")]:
        sel=data[(data.design==design)&(data.endpoint==ep)].set_index("panel")
        ordered.extend([sel.loc[k] for k in (joint,"RNA_only","protein_only")])
    for y,row in zip(ys,ordered):
        ax.barh(y,31,color="#f3f2ee",height=.62)
        ax.barh(y,row.testable,color=pale,height=.62)
        ax.barh(y,row.confirmed,color=ink,height=.62)
        ax.text(32,y,f"{int(row.confirmed)} / {int(row.testable)}",va="center",fontsize=8.5)
    ax.text(0,6.47,"Prespecified panels",fontsize=8.6,weight="bold")
    ax.text(0,2.66,"Exploratory panels",fontsize=8.6,weight="bold")
    ax.set(yticks=ys,yticklabels=labels,xticks=[0,10,20,31],xlim=(0,41),ylim=(-.65,7.05))
    ax.set_title(title,loc="left",weight="bold")
    ax.set_xlabel("Nominations (31 slots per row)")
    ax.tick_params(axis="y",length=0)
    ax.spines[["top","right","left"]].set_visible(False)

ax=fig.add_subplot(grid[1,:])
for i,(ep,model,design) in enumerate([("RNA","RNA_only","primary"),("RNA","RNA_only","exploratory"),("protein","protein_only","primary"),("protein","protein_only","exploratory")]):
    r=diff[(diff.endpoint==ep)&(diff.comparator==model)&(diff.design==design)].iloc[0]
    ax.plot([r.stability_low,r.stability_high],[3-i,3-i],color=ink,lw=1.8)
    ax.scatter(r.difference,3-i,color=orange,marker="o",s=32,zorder=5)
ax.axvline(0,ls="--",color="0.45",lw=.9)
ax.set(yticks=range(3,-1,-1),yticklabels=["RNA · joint only","RNA · joint overall","Protein · joint only","Protein · joint overall"],
       ylim=(-.6,3.6),xlim=(-20,5),xticks=[-20,-15,-10,-5,0,5])
ax.set_xlabel("Additional confirmations versus selection by the matching assay")
ax.set_title("C   Paired differences and 95% patient-resampling ranges",loc="left",weight="bold")
ax.tick_params(axis="y",length=0)
ax.spines[["left","right","top"]].set_visible(False)
fig.legend(handles=[Patch(facecolor=ink,label="Confirmed"),Patch(facecolor=pale,label="Testable, unconfirmed"),Patch(facecolor="#f3f2ee",edgecolor="#c9c8c2",label="Untestable")],
           loc="outside lower center",ncol=3,frameon=False)
for ext in ("svg","png"):
    fig.savefig(root/("glioma-followup."+ext),dpi=300,metadata={"Date":None} if ext=="svg" else None)
plt.close(fig)
