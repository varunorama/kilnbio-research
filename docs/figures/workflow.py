"""Draw a capability-only KilnBio workflow with ReportLab; no study data required."""
from pathlib import Path
import base64, csv, json, math
from xml.sax.saxutils import escape
from reportlab.pdfgen.canvas import Canvas
from reportlab.lib.colors import HexColor
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
import reportlab

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
FONT=Path(reportlab.__file__).parent/"fonts"
for face,file in [("Vera","Vera.ttf"),("Vera-Bold","VeraBd.ttf")]:
    pdfmetrics.registerFont(TTFont(face,str(FONT/file)))
ink="#123237"; teal="#006A68"; orange="#A14725"; muted="#435D60"
line="#8BA8A6"; grid="#C8D7D5"; pale="#ECF5F3"; warm="#FAF1E8"; white="#FFFFFF"
paper="#FFFFFF"; blue="#315D77"

class Graphic:
    def __init__(self,name,height,summary):
        self.h=height; self.name=name;self.scale=493.23/700
        self.c=Canvas(str(OUT/(name+".pdf")),pagesize=(700*self.scale,height*self.scale),invariant=1)
        self.c.scale(self.scale,self.scale); self.c.setAuthor("");self.c.setTitle(name.replace("_"," "))
        self.body=[f'<svg xmlns="http://www.w3.org/2000/svg" width="174mm" height="{174*height/700:.2f}mm" viewBox="0 0 700 {height}" role="img" aria-labelledby="title desc"><title id="title">{escape(name.replace("_"," "))}</title><desc id="desc">{escape(summary)}</desc>']
        if name=="Fig1":
            faces=[("normal","Vera.ttf"),("bold","VeraBd.ttf")]
            styles="".join("@font-face{font-family:KilnVera;font-weight:"+weight+";src:url(data:font/ttf;base64,"+base64.b64encode((FONT/fn).read_bytes()).decode("ascii")+") format('truetype');}" for weight,fn in faces)
            self.body.append("<defs><style>"+styles+"</style></defs>")
        self.rect(0,0,700,height,paper)
    def rect(self,x,y,w,h,fill,edge=None,r=0):
        c=self.c;c.setFillColor(HexColor(fill))
        if edge:c.setStrokeColor(HexColor(edge));c.setLineWidth(1)
        c.roundRect(x,self.h-y-h,w,h,r,stroke=bool(edge),fill=1)
        self.body.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="{fill}" stroke="{edge or "none"}"/>')
    def text(self,x,y,s,z=14,bold=False,col=ink,align="left",maxw=None):
        s=str(s); f="Vera-Bold" if bold else "Vera"
        a=pdfmetrics.stringWidth(s,f,z)
        if maxw and a>maxw+0.01:raise ValueError(f"Overlong {s!r}: {a:.1f} > {maxw}")
        left=x-(a if align=="right" else a/2 if align=="center" else 0)
        if left < -0.1 or left+a>700.1:raise ValueError(f"Outside artboard: {s}")
        c=self.c;c.setFont(f,z);c.setFillColor(HexColor(col));c.drawString(left,self.h-y,s)
        self.body.append(f'<text x="{left:.3f}" y="{y}" font-family="KilnVera,DejaVu Sans,Arial,sans-serif" font-weight="{ "bold" if bold else "normal" }" font-size="{z}" fill="{col}">{escape(s)}</text>')
    def path(self,pts,col=line,width=1.6,dash=None,closed=False,fill=None):
        c=self.c;c.saveState();c.setStrokeColor(HexColor(col));c.setLineWidth(width);c.setLineCap(1);c.setLineJoin(1)
        if dash:c.setDash(dash)
        if fill:c.setFillColor(HexColor(fill))
        p=c.beginPath();p.moveTo(pts[0][0],self.h-pts[0][1])
        for x,y in pts[1:]:p.lineTo(x,self.h-y)
        if closed:p.close()
        c.drawPath(p,stroke=1,fill=bool(fill));c.restoreState()
        ds=f' stroke-dasharray="{dash[0]} {dash[1]}"' if dash else ""
        d=" ".join(("M" if i==0 else "L")+f"{x:.3f} {y:.3f}" for i,(x,y) in enumerate(pts))+(" Z" if closed else "")
        self.body.append(f'<path d="{d}" fill="{fill or "none"}" stroke="{col}" stroke-width="{width}" stroke-linecap="round" stroke-linejoin="round"{ds}/>')
    def dot(self,x,y,r,fill,edge=None,width=1.4):
        c=self.c;c.setFillColor(HexColor(fill));c.setStrokeColor(HexColor(edge or fill));c.setLineWidth(width)
        c.circle(x,self.h-y,r,stroke=bool(edge),fill=1)
        self.body.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}" stroke="{edge or "none"}" stroke-width="{width}"/>')
    def arrow(self,pts,col=teal,width=2,head=5):
        self.path(pts,col,width)
        x,y=pts[-1];xx,yy=pts[-2];ang=math.atan2(y-yy,x-xx)
        one=(x-head*math.cos(ang-.52),y-head*math.sin(ang-.52));two=(x-head*math.cos(ang+.52),y-head*math.sin(ang+.52))
        self.path([one,(x,y),two],col,width)
    def close(self,svg=False):
        self.c.showPage();self.c.save()
        if svg:(OUT/(self.name+".svg")).write_text("\n".join(self.body+["</svg>"]),encoding="utf-8")

def assay_glyph(g,x,y):
    # Aligned segments, an accessibility peak and per-site marks depict three assay representations.
    for start,end,dy in [(0,14,0),(17,25,0),(5,19,4)]:g.path([(x+start,y+dy),(x+end,y+dy)],teal,2)
    g.path([(x,y+17),(x+6,y+17),(x+10,y+12),(x+15,y+17),(x+19,y+9),(x+22,y+17),(x+26,y+17)],blue,1.65)
    g.path([(x,y+25),(x+26,y+25)],line,1.2)
    for i in [3,11,21]:g.dot(x+i,y+25,2.4,orange)
def match_glyph(g,x,y):
    for dy in (0,11,22):
        g.dot(x,y+dy,2.8,teal);g.dot(x+20,y+dy,2.8,blue)
    for dy in (0,11):g.path([(x+4,y+dy),(x+16,y+dy)],line,1.5)
    g.path([(x+5,y+22),(x+15,y+22)],line,1.2,[2,2])
def branch_glyph(g,x,y):
    g.path([(x,y+11),(x+9,y+11),(x+9,y),(x+22,y)],teal,1.8)
    g.path([(x+9,y+11),(x+9,y+22),(x+22,y+22)],teal,1.8)
    for xx,yy in [(x,y+11),(x+22,y),(x+22,y+22)]:g.dot(xx,yy,3,orange if yy==y else teal)
def mini(g,x,y,kind):
    if kind=="reads":
        for i,a,b in [(0,0,13),(1,6,23),(2,1,9)]:g.path([(x+a,y+5*i),(x+b,y+5*i)],teal,2.4)
    elif kind=="cells":
        for a,b in [(3,4),(13,2),(9,12),(19,11)]:g.dot(x+a,y+b,3.1,white,teal,1.4)
    elif kind=="peaks":
        g.path([(x,y+13),(x+4,y+12),(x+8,y+3),(x+12,y+11),(x+16,y+7),(x+20,y),(x+24,y+13)],blue,1.7)
    elif kind=="sites":
        g.path([(x,y+11),(x+24,y+11)],teal,1.1)
        for i in (2,10,19):
            g.path([(x+i,y+11),(x+i,y+5)],teal,1.2);g.dot(x+i,y+4,2.5,white if i==10 else teal,teal,1)
    elif kind=="matrix":
        for i in range(4):
            for j in range(3):g.rect(x+i*5.5,y+5*j,3.6,3.3,[blue,line,teal][(i+j)%3])
    elif kind=="amount":
        for i,h in enumerate([5,12,8,15]):g.rect(x+i*6,y+14-h,3.7,h,teal)
    elif kind=="mods":
        for j in range(2):
            g.path([(x,y+3+7*j),(x+24,y+3+7*j)],line,1)
            for i in range(3):g.dot(x+3+i*9,y+3+7*j,2.25,white if (i+j)%2 else teal,teal,1)
    else:
        for i in range(4):
            g.rect(x+i*6,y+3,4.3,8,orange if i==2 else grid,r=1)
        g.path([(x+14,y),(x+14,y-2)],orange,1.6)

def draw(branded=True):
    g=Graphic("Fig1",740,"KilnBio integrates assay-specific measurements through matched biological units and explicit gene or region maps. Optional genomic and protein annotations remain distinct from measurements. Outputs support research comparisons and experimental follow-up, not clinical diagnosis.")
    if branded:
        icon=ROOT/"docs/brand/kilnbio-icon.png"
        g.c.drawImage(str(icon),18,g.h-58,50,50,mask="auto")
        raw=base64.b64encode(icon.read_bytes()).decode()
        g.body.append(f'<image x="18" y="8" width="50" height="50" href="data:image/png;base64,{raw}"/>')
    x=79 if branded else 18
    g.text(x,34,"KilnBio builds context across omics",23,True,maxw=600)
    g.text(x,54,"Research use  /  Process, connect, compare and test molecular evidence",12.7,False,muted,maxw=603)
    def head(y,num,title):
        g.dot(27,y-5,9,teal);g.text(27,y-1.2,num,10,True,white,"center")
        g.text(43,y,title,13.2,True,teal)
    head(84,"1","Start with raw reads or existing molecular measurements")
    entries=[(18,"Sequence reads",[
        ("Bulk RNA-seq","STAR + RSEM or Salmon; gene counts","reads"),
        ("10x single-cell RNA-seq","STARsolo, cell QC and donor pseudobulk","cells"),
        ("ATAC-seq and ChIP-seq","Peak callers and shared-region counts","peaks"),
        ("WGBS and RRBS","Bisulfite alignment and cytosine counts","sites"),
        ("Autosomal exome capture","Sarek; selectable germline callers","variants")]),
        (358,"Existing measurements",[
        ("RNA, peaks and per-feature effects","Count matrices or reviewed statistics","matrix"),
        ("Protein and metabolite abundance","Olink NPX or processed assay matrices","amount"),
        ("DNA modifications and arrays","ONT modBAM, biomodal or beta values","mods"),
        ("Genomic states and CH candidates","Called variants and explicit callability","variants"),
        ("Custom and commercial panels","Exome or methylome target BEDs","sites")])]
    for x,title,rows in entries:
        g.rect(x,95,324,203,pale if x==18 else warm,grid,7)
        (assay_glyph if x==18 else match_glyph)(g,x+13,105)
        g.text(x+52,119,title,17,True,maxw=265)
        for j,(a,b,kind) in enumerate(rows):
            y=148+j*31;mini(g,x+13,y-6,kind)
            g.text(x+45,y,a,12.3,True,maxw=269)
            g.text(x+45,y+14,b,11.0,False,muted,maxw=269)
    g.path([(180,298),(180,308),(520,308),(520,298)],teal,1.8)
    g.arrow([(350,308),(666,308),(666,336)],width=1.8)
    head(330,"2","Connect measurements under one research question")
    g.rect(18,340,664,76,pale,grid,7)
    for x,title,lines in [
        (32,"Match biological units",["Donor, specimen, tissue and time;","track missing assays and coverage"]),
        (257,"Link molecular features",["Alleles, CpGs, regions and genes;","record mapping assumptions"]),
        (465,"Model each assay",["Preserve units and uncertainty;","fit a shared contrast where valid"])]:
        g.text(x,361,title,12.8,True,maxw=204)
        for j,t in enumerate(lines):g.text(x,381+j*16,t,11.2,False,muted,maxw=207)
    for x in [241,450]:g.path([(x,352),(x,403)],grid,1)
    g.arrow([(666,416),(666,444)],width=1.7)
    head(439,"3","Add optional context while keeping evidence types distinct")
    g.rect(18,448,664,64,warm,grid,7)
    g.text(32,468,"VEP consequences  |  ClinVar assertions  |  AlphaMissense scores",12.7,True,maxw=636)
    g.text(32,486,"AlphaFold Database structures  |  AlphaGenome predictions",12.7,True,maxw=636)
    g.text(32,503,"Exact allele / transcript / protein mapping; annotations do not replace measured effects.",11.3,False,muted,maxw=636)
    g.path([(666,512),(666,552),(125,552)],teal,1.5)
    for x in [125,350,571]:g.arrow([(x,552),(x,561)],width=1.6,head=4)
    head(547,"4","Ask what the integrated evidence supports")
    for x,title,rt,body,end in [
        (18,"Compare analysis choices","Same specimens and features",["Differential methods, peak callers","and assay combinations"],"Agreement and sensitivity"),
        (243,"Prioritize hypotheses","Gene or regulatory context",["Convergent or conflicting effects;","explicit missing evidence"],"Candidates for follow-up"),
        (468,"Test generalization","Held-out cohorts or materials",["Replication, prediction and","preservation of biological contrasts"],"Yield, error and exclusions")]:
        g.rect(x,562,214,91,white,grid,6);g.rect(x+1,563,212,4,teal if x!=468 else orange,r=1)
        g.text(x+10,584,title,12.6,True,maxw=194)
        g.text(x+10,600,rt,10.9,True,teal,maxw=195)
        for j,t in enumerate(body):g.text(x+10,616+13*j,t,10.65,False,muted,maxw=195)
        g.text(x+10,645,end,11.1,True,maxw=195)
    g.path([(18,666),(682,666)],grid,1)
    g.text(19,685,"RESEARCH EXAMPLE",10.3,True,teal)
    g.text(169,685,"Regulatory context in K562 cells",13.4,True)
    for x,a,b in [(19,"Accessibility + H3K27ac","+ measured CpG methylation"),
                  (254,"Map to enhancer–gene pairs","and compare eligible rankings"),
                  (487,"Freeze candidate links","then test CRISPR support")]:
        g.text(x,709,a,11.7,True,maxw=200);g.text(x,726,b,11.3,False,muted,maxw=202)
    for x in (219,450):g.arrow([(x,714),(x+22,714)],width=1.6,head=4)
    g.close(svg=True)

if __name__=="__main__":
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir",type=Path,required=True);p.add_argument("--no-icon",action="store_true")
    a=p.parse_args();OUT=a.output_dir.resolve();OUT.mkdir(parents=True,exist_ok=True);draw(not a.no_icon)
