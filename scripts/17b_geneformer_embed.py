#!/usr/bin/env python
# 17b — Geneformer zero-shot 嵌入 grade 分层肿瘤细胞，检验基础模型是否捕获侵袭性轴
import os, scipy.io, numpy as np, pandas as pd, anndata as ad
ED="results/scrna/gf_export"; OUT="results/scrna/gf_out"; os.makedirs(OUT,exist_ok=True)

# 1) 组 AnnData（cells x genes）
M=scipy.io.mmread(os.path.join(ED,"counts.mtx")).tocsr()   # genes x cells
genes=pd.read_csv(os.path.join(ED,"genes.csv")); cells=pd.read_csv(os.path.join(ED,"cells.csv"))
X=M.T.tocsr()
adata=ad.AnnData(X=X, obs=cells.set_index("cell"), var=genes.set_index("symbol"))
# Preserve stable cell- and patient-level identifiers through Geneformer
# tokenisation. EmbExtractor may reorder rows, so downstream analyses must not
# assume that the embedding CSV has the same row order as cells.csv.
adata.obs["cell_id"]=adata.obs.index.astype(str)
adata.var["ensembl_id"]=genes["ensembl_id"].values
adata.obs["n_counts"]=np.asarray(X.sum(1)).ravel()
adata.obs["grade"]=adata.obs["grade"].astype(str)
h5=os.path.join(ED,"tumor.h5ad"); adata.write_h5ad(h5)
print("anndata:",adata.shape,"-> ",h5)

# 2) Tokenize (Geneformer V2: special_token=True, input 2048)
from geneformer import TranscriptomeTokenizer
tk=TranscriptomeTokenizer({"grade":"grade","AggrScore":"AggrScore",
                            "gsm":"gsm","cell_id":"cell_id"}, nproc=2,
                          model_input_size=2048, special_token=True)
tk.tokenize_data(ED, OUT, "tumor", file_format="h5ad")
print("tokenized ->", os.path.join(OUT,"tumor.dataset"))

# 3) 嵌入提取（pretrained，zero-shot）
from geneformer import EmbExtractor
mdir="results/scrna/gf_model_V2-104M"
emb=EmbExtractor(model_type="Pretrained", num_classes=0, emb_mode="cell",
                 emb_layer=-1, max_ncells=None, forward_batch_size=8, nproc=2,
                 emb_label=["grade","AggrScore","gsm","cell_id"])
df=emb.extract_embs(mdir, os.path.join(OUT,"tumor.dataset"), OUT, "tumor_emb")
print("emb df:",df.shape)

# 4) 分析：嵌入是否按 grade 组织 + 与程序评分相关
embcols=[c for c in df.columns if isinstance(c,int) or str(c).isdigit()]
E=df[embcols].values.astype(float)
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
Ez=StandardScaler().fit_transform(E)
pca=PCA(n_components=10,random_state=42).fit(Ez); P=pca.transform(Ez)
gmap={"I":1,"II":2,"III":3}; g=df["grade"].map(gmap).values
from scipy.stats import spearmanr
# 与 grade 相关性最强的 PC
rho_g=[spearmanr(P[:,i],g).correlation for i in range(P.shape[1])]
best=int(np.argmax(np.abs(rho_g)))
sil=silhouette_score(Ez, df["grade"].values)
res=pd.DataFrame({"PC":range(1,11),"spearman_vs_grade":np.round(rho_g,3)})
res.to_csv(os.path.join(OUT,"geneformer_emb_grade_corr.csv"),index=False)
with open(os.path.join(OUT,"geneformer_summary.txt"),"w") as f:
    f.write(f"cells={E.shape[0]} embdim={E.shape[1]}\n")
    f.write(f"silhouette(grade)={sil:.3f}\n")
    f.write(f"best PC vs grade: PC{best+1} rho={rho_g[best]:.3f}\n")
    if "AggrScore" in df: f.write(f"PC{best+1} vs AggrScore rho={spearmanr(P[:,best],pd.to_numeric(df['AggrScore'],errors='coerce')).correlation:.3f}\n")
print("silhouette(grade)=",round(sil,3),"bestPC vs grade rho=",round(rho_g[best],3))

# 5) UMAP of embedding colored by grade + score
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    import umap
    U=umap.UMAP(random_state=42).fit_transform(Ez)
    fig,axs=plt.subplots(1,2,figsize=(11,4.6))
    for gi,gc in zip(["I","II","III"],["#3b6ea5","#caa54b","#a14b3d"]):
        m=df["grade"].values==gi; axs[0].scatter(U[m,0],U[m,1],s=5,c=gc,label=gi,alpha=.6)
    axs[0].legend(title="WHO grade"); axs[0].set_title(f"Geneformer zero-shot embedding\nsilhouette(grade)={sil:.2f}")
    sc=axs[1].scatter(U[:,0],U[:,1],s=5,c=pd.to_numeric(df["AggrScore"],errors="coerce"),cmap="RdBu_r")
    plt.colorbar(sc,ax=axs[1]); axs[1].set_title("colored by bulk aggressiveness score")
    plt.tight_layout(); plt.savefig(os.path.join("results/figures_pub","fig_geneformer_embedding.png"),dpi=150)
    plt.savefig(os.path.join("results/figures_pub","fig_geneformer_embedding.pdf"))
    print("saved fig_geneformer_embedding")
except Exception as e:
    print("umap/plot skipped:",str(e)[:120])
print("EMBED-DONE")
