#!/usr/bin/env python3
"""Generate the meningioma graphical abstract as a 4:3 vector SVG.
Content truth source: docs/GA_design_brief.md ; mirrors locked raster docs/GA/round2.png.
All numbers come from the brief whitelist; wording obeys the causal-language red lines."""
import random, math
from pathlib import Path
random.seed(42)
W, H = 1200, 900

# Okabe-Ito colourblind-safe palette
BLUE="#0072B2"; ORANGE="#E69F00"; VERM="#D55E00"; GREEN="#009E73"
PURPLE="#CC79A7"; SKY="#56B4E9"; GREY="#9C9C9C"; DARK="#222222"; LIGHT="#E9E9E9"

S=[]
def esc(t): return t.replace("&","&amp;").replace("<","&lt;").replace(">","&gt;")
def rect(x,y,w,h,rx=0,fill="none",stroke="none",sw=1,op=1,dash=None):
    d=f' stroke-dasharray="{dash}"' if dash else ''
    S.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" opacity="{op}"{d}/>')
def circ(cx,cy,r,fill="none",stroke="none",sw=1,op=1):
    S.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" opacity="{op}"/>')
def line(x1,y1,x2,y2,stroke=DARK,sw=1,dash=None,cap="round"):
    d=f' stroke-dasharray="{dash}"' if dash else ''
    S.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{stroke}" stroke-width="{sw}" stroke-linecap="{cap}"{d}/>')
def path(d,fill="none",stroke="none",sw=1):
    S.append(f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" stroke-linejoin="round" stroke-linecap="round"/>')
def text(x,y,t,size=14,fill=DARK,weight="normal",anchor="start",style="normal",ff="Arial"):
    S.append(f'<text x="{x}" y="{y}" font-family="{ff}" font-size="{size}" font-weight="{weight}" font-style="{style}" fill="{fill}" text-anchor="{anchor}">{esc(t)}</text>')

def arrow(x1,y1,x2,y2,stroke=GREY,sw=3):
    line(x1,y1,x2,y2,stroke=stroke,sw=sw)
    ang=math.atan2(y2-y1,x2-x1); L=10
    for s in (0.5,-0.5):
        ax=x2-L*math.cos(ang-s); ay=y2-L*math.sin(ang-s)
        line(x2,y2,ax,ay,stroke=stroke,sw=sw)

def badge(x,y,n,color):
    rect(x,y,30,30,rx=7,fill=color)
    text(x+15,y+21,str(n),size=18,fill="#ffffff",weight="bold",anchor="middle")

def panel(x,y,w,h,color):
    rect(x,y,w,h,rx=14,fill="#ffffff",stroke=color,sw=2.5)

# ---------- background ----------
rect(0,0,W,H,fill="#ffffff")
# ---------- title ----------
text(W/2,46,"Human meningioma grade-associated aggressiveness program",size=30,weight="bold",anchor="middle")

PX1,PX2=20,610; PW=570
PY1,PY2=70,448; PH=370

# ================= PANEL 1 : Bulk program definition =================
panel(PX1,PY1,PW,PH,BLUE); badge(PX1+16,PY1+16,1,BLUE)
text(PX1+56,PY1+38,"Bulk program definition",size=20,weight="bold")
rect(PX1+395,PY1+18,150,26,rx=13,fill="#ffffff",stroke=BLUE,sw=1.5)
text(PX1+470,PY1+36,"human meningioma",size=12,fill=BLUE,anchor="middle")
# three dataset cards
cards=[(BLUE,"GSE136661 (n≈160)","discovery"),(SKY,"GSE16581 n=68","validation"),(GREEN,"GSE74385 n=62","validation")]
cx0=PX1+24; cw=160; gap=14
for i,(col,l1,l2) in enumerate(cards):
    cx=cx0+i*(cw+gap)
    rect(cx,PY1+62,cw,52,rx=8,fill="#ffffff",stroke=col,sw=1.8)
    rect(cx+12,PY1+76,24,24,rx=3,fill=col)
    rect(cx+46,PY1+80,100,6,rx=3,fill=LIGHT); rect(cx+46,PY1+94,80,6,rx=3,fill=LIGHT)
    text(cx+cw/2,PY1+132,l1,size=12,anchor="middle",weight="bold")
    text(cx+cw/2,PY1+148,l2,size=12,anchor="middle",fill=GREY)
# connector down
midx=PX1+24+cw/2 + (cw+gap)  # under middle card
line(cx0+cw/2,PY1+158,cx0+2*(cw+gap)+cw/2,PY1+158,stroke=GREY,sw=2)
arrow(midx, PY1+158, midx, PY1+184, stroke=GREY, sw=2.5)
# program box
bx,by,bw,bh=PX1+70,PY1+188,360,90
rect(bx,by,bw,bh,rx=10,fill="#F2FBF8",stroke=GREEN,sw=2)
text(bx+bw/2,by+30,"200-gene program",size=20,weight="bold",fill=GREEN,anchor="middle")
text(bx+bw/2,by+52,"up  OXPHOS / ribosome",size=13,anchor="middle")
text(bx+bw/2,by+72,"down  synaptic / Wnt",size=13,anchor="middle")
# small bar icon to the right
ix=bx+bw+18
for j,hh in enumerate((16,26,36)):
    rect(ix+j*14,by+58-hh,9,hh,rx=2,fill=GREEN)
circ(ix+21,by+45,42,fill="none",stroke=GREEN,sw=2)

# ================= PANEL 2 : Grade & recurrence =================
panel(PX2,PY1,PW,PH,ORANGE); badge(PX2+16,PY1+16,2,ORANGE)
text(PX2+56,PY1+38,"Grade & recurrence association",size=20,weight="bold")
# boxplot axes
ax,ay=PX2+40,PY1+78; aw,ah=250,210
line(ax,ay,ax,ay+ah,stroke=DARK,sw=2)      # y axis
line(ax,ay+ah,ax+aw,ay+ah,stroke=DARK,sw=2) # x axis
S.append(f'<text x="{ax-12}" y="{ay+ah/2}" font-family="Arial" font-size="12" fill="{DARK}" text-anchor="middle" transform="rotate(-90 {ax-12} {ay+ah/2})">program score</text>')
# remove the non-rotated dup (overwrite by drawing white)—simpler: ignore, the rotated one drawn after sits on top
boxes=[(BLUE,"WHO I",150,40),(ORANGE,"WHO II",95,55),(VERM,"WHO III",55,70)]
bw2=46
for i,(col,lab,top,boxh) in enumerate(boxes):
    bxc=ax+45+i*70
    ytop=ay+top
    # whiskers
    line(bxc,ytop-22,bxc,ytop+boxh+22,stroke=DARK,sw=1.5)
    line(bxc-12,ytop-22,bxc+12,ytop-22,stroke=DARK,sw=1.5)
    line(bxc-12,ytop+boxh+22,bxc+12,ytop+boxh+22,stroke=DARK,sw=1.5)
    rect(bxc-bw2/2,ytop,bw2,boxh,rx=2,fill=col,op=0.9,stroke=DARK,sw=1)
    line(bxc-bw2/2,ytop+boxh*0.45,bxc+bw2/2,ytop+boxh*0.45,stroke="#ffffff",sw=1.5)
    text(bxc,ay+ah+18,lab,size=12,anchor="middle")
# rising arrow
arrow(ax+40,ay+30,ax+aw-10,ay-6,stroke=VERM,sw=2.5)
text(ax+50,ay-12,"program score rises",size=12,fill=VERM,style="italic")
# text box
tx,ty,tw,th=PX2+310,PY1+90,235,205
rect(tx,ty,tw,th,rx=10,fill="#FFF8EC",stroke=ORANGE,sw=2)
text(tx+tw/2,ty+26,"associated with recurrence",size=14,weight="bold",fill=VERM,anchor="middle")
lines2=["GSE16581 grade ANOVA p=1.7e-6","GSE16581 recurrence p=4.8e-6","","adjusted, retrospective:","GSE16581 OR/1-SD 14.74 (2.77–224.07)","GSE74385 recurrence-only 5.18 (1.40–25.83)"]
yy=ty+52
for ln in lines2:
    if ln: text(tx+14,yy,ln,size=12.5)
    yy+=24

# ================= PANEL 3 : Single-cell basis =================
panel(PX1,PY2,PW,PH,GREEN); badge(PX1+16,PY2+16,3,GREEN)
text(PX1+56,PY2+38,"Single-cell cellular basis",size=20,weight="bold")
text(PX1+24,PY2+64,"GSE183655",size=13,weight="bold",fill=GREEN)
text(PX1+24,PY2+82,"GSE206647  163,897 cells",size=13,weight="bold",fill=GREEN)
text(PX1+24,PY2+100,"16 grade-labelled patients",size=12,fill=DARK)
# UMAP scatter clusters
clusters=[(BLUE,90,200),(SKY,150,168),(GREEN,140,250),(ORANGE,210,165),(PURPLE,225,235),(VERM,95,255),(GREY,80,300)]
ox,oy=PX1+20,PY2+110
for col,cxr,cyr in clusters:
    for _ in range(26):
        a=random.random()*6.28; r=random.random()**0.5*22
        circ(ox+cxr+r*math.cos(a),oy+cyr*0.55+r*math.sin(a),2.4,fill=col,op=0.85)
# highlight tumour cluster (green)
circ(ox+140,oy+250*0.55,30,fill="none",stroke=GREEN,sw=2.2)
text(ox+140,oy+250*0.55+52,"tumour/tumour-like cells",size=12,weight="bold",fill=GREEN,anchor="middle")
# arrow to line plot
arrow(PX1+260,PY2+200,PX1+300,PY2+200,stroke=GREY,sw=2.5)
# line plot per-patient trajectories
lx,ly=PX1+330,PY2+120; lw,lh=200,170
line(lx,ly,lx,ly+lh,stroke=DARK,sw=2); line(lx,ly+lh,lx+lw,ly+lh,stroke=DARK,sw=2)
S.append(f'<text x="{lx-10}" y="{ly+lh/2}" font-family="Arial" font-size="12" fill="{DARK}" text-anchor="middle" transform="rotate(-90 {lx-10} {ly+lh/2})">program score</text>')
xs=[lx+30,lx+100,lx+170]
for k in range(6):  # grey patient lines
    base=[random.uniform(0.15,0.45),random.uniform(0.45,0.7),random.uniform(0.7,0.95)]
    pts=[(xs[i],ly+lh-base[i]*lh) for i in range(3)]
    path(f'M {pts[0][0]} {pts[0][1]} L {pts[1][0]} {pts[1][1]} L {pts[2][0]} {pts[2][1]}',stroke=GREY,sw=1.3)
    for px,py in pts: circ(px,py,2.5,fill="#ffffff",stroke=GREY,sw=1.2)
mean=[(xs[0],ly+lh-0.25*lh),(xs[1],ly+lh-0.6*lh),(xs[2],ly+lh-0.85*lh)]
path(f'M {mean[0][0]} {mean[0][1]} L {mean[1][0]} {mean[1][1]} L {mean[2][0]} {mean[2][1]}',stroke=GREEN,sw=3.5)
for px,py in mean: circ(px,py,5,fill=GREEN)
for i,lab in enumerate(("WHO I","WHO II","WHO III")):
    text(xs[i],ly+lh+18,lab,size=11,anchor="middle")
text(lx+lw/2,PY2+PH-16,"pseudobulk rho=0.87, p=1.2e-5",size=13,weight="bold",fill=GREEN,anchor="middle")

# ================= PANEL 4 : Foundation-model representation =================
panel(PX2,PY2,PW,PH,PURPLE); badge(PX2+16,PY2+16,4,PURPLE)
text(PX2+56,PY2+38,"Foundation-model representation",size=19,weight="bold")
rect(PX2+360,PY2+18,185,26,rx=13,fill="#ffffff",stroke=PURPLE,sw=1.5)
text(PX2+452,PY2+36,"Geneformer V2-104M zero-shot",size=10.5,fill=PURPLE,anchor="middle")
# left subpanel: embedding scatter
ex,ey,ew,eh=PX2+18,PY2+58,300,230
rect(ex,ey,ew,eh,rx=10,fill="#FBF4F8",stroke=PURPLE,sw=2)
text(ex+ew/2,ey+24,"patient-level representation",size=13,weight="bold",fill=PURPLE,anchor="middle")
grp=[(BLUE,70,150,"WHO I"),(ORANGE,150,110,"WHO II"),(VERM,225,80,"WHO III")]
for col,cxr,cyr,lab in grp:
    for _ in range(34):
        a=random.random()*6.28; r=random.random()**0.5*26
        circ(ex+cxr+r*math.cos(a),ey+cyr+r*math.sin(a)*0.8,2.5,fill=col,op=0.8)
arrow(ex+45,ey+eh-40,ex+ew-40,ey+50,stroke=PURPLE,sw=2.5)
for col,cxr,cyr,lab in grp:
    text(ex+cxr,ey+eh-12,lab,size=10.5,fill=col,anchor="middle",weight="bold")
text(ex+ew/2,ey+eh+0,"",size=1)
# rho text under left subpanel
text(ex+ew/2,PY2+PH-16,"primary |rho| 0.684/0.665   perm p=0.0516",size=12.5,weight="bold",fill=PURPLE,anchor="middle")
# right subpanel: honest benchmark bars
hx,hy,hw,hh=PX2+332,PY2+58,212,200
rect(hx,hy,hw,hh,rx=10,fill="#F2F7FC",stroke=BLUE,sw=2)
text(hx+hw/2,hy+24,"honest side-by-side benchmark",size=11.5,weight="bold",fill=BLUE,anchor="middle")
base_y=hy+150
# bars scaled: AUC mapped from 0.5..1.0 to height
def bh(auc): return (auc-0.5)/0.5*95
b1x=hx+45; b2x=hx+125; bwid=46
rect(b1x,base_y-bh(0.99),bwid,bh(0.99),rx=3,fill=BLUE)
rect(b2x,base_y-bh(0.977),bwid,bh(0.977),rx=3,fill=PURPLE)
line(hx+20,base_y,hx+hw-20,base_y,stroke=DARK,sw=1.5)
text(b1x+bwid/2,base_y+18,"classical",size=11,anchor="middle",fill=BLUE,weight="bold")
text(b1x+bwid/2,base_y+33,"AUC 0.99",size=11,anchor="middle",fill=BLUE)
text(b2x+bwid/2,base_y+18,"FM",size=11,anchor="middle",fill=PURPLE,weight="bold")
text(b2x+bwid/2,base_y+33,"AUC 0.98",size=11,anchor="middle",fill=PURPLE)
text(hx+hw/2,hy+hh-6,"n=16; descriptive",size=11.5,anchor="middle",style="italic")

# ================= take-home band =================
rect(20,824,W-40,62,rx=12,fill="#F6F6F6",stroke=GREY,sw=1.5)
th1="A reproducible cross-platform grade-associated aggressiveness program is associated with meningioma recurrence,"
th2="resolved in tumour/tumour-like cells, with an aggregation-dependent Geneformer benchmark."
text(W/2,852,th1,size=15,weight="bold",anchor="middle")
text(W/2,873,th2,size=15,weight="bold",anchor="middle")

svg=f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">\n'+"\n".join(S)+"\n</svg>\n"
root = Path(__file__).resolve().parents[1]
(root / "docs" / "GA" / "GA_final.svg").write_text(svg, encoding="utf-8")
print("wrote GA_final.svg", len(svg),"bytes")
