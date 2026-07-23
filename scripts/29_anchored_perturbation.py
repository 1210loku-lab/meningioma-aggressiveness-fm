#!/usr/bin/env python
# 29 — 锚定 in-silico perturbation：敲除哪些 top 程序基因把肿瘤细胞从"高侵袭态"推向"低侵袭态"
# 正规 Geneformer goal-state-shift 工作流 + 对标经典 RF 重要性
# macOS: 强制 fork 启动方式，绕过 datasets 5.0 + spawn 的 multiprocess Manager EOFError（须在 import datasets 前）
import multiprocessing as _mp
try: _mp.set_start_method("fork", force=True)
except RuntimeError: pass
try:
    import multiprocess as _mpd; _mpd.set_start_method("fork", force=True)
except Exception: pass
import os, pickle, numpy as np, pandas as pd
from datasets import load_from_disk
from geneformer import EmbExtractor, InSilicoPerturber, InSilicoPerturberStats
np.random.seed(42)
OUT="results/scrna/gf_anchor"; os.makedirs(OUT,exist_ok=True)
MODEL="results/scrna/gf_model_V2-104M"
log=open("results/deg/P2_anchored_perturbation.txt","w")
def w(s): print(s); log.write(str(s)+"\n"); log.flush()

# 1) tokenized dataset + high/low 侵袭态（复用已存则跳过）
ds_path=os.path.join(OUT,"tumor_states.dataset")
import geneformer, collections
if os.path.isdir(ds_path):
    d=load_from_disk(ds_path); w("reuse existing tumor_states.dataset")
else:
    d=load_from_disk("results/scrna/gf_out/tumor.dataset")
    ag=np.array(d["AggrScore"]); q1,q2=np.quantile(ag,[1/3,2/3])
    state=np.where(ag>=q2,"high",np.where(ag<=q1,"low","mid"))
    d=d.add_column("aggr_state", state.tolist()); d=d.filter(lambda e: e["aggr_state"]!="mid")
    d.save_to_disk(ds_path)
w(f"high/low cells: {sum(np.array(d['aggr_state'])=='high')}/{sum(np.array(d['aggr_state'])=='low')}")

# 2) genes_to_perturb = 程序基因中在细胞里实际被 token 化、出现频率最高的 top30
tokd=pickle.load(open(os.path.join(os.path.dirname(geneformer.__file__),"token_dictionary_gc104M.pkl"),"rb"))
drv=pd.read_csv("results/deg/classical_drivers_ranked.csv")
freq=collections.Counter()
for ids in d["input_ids"]: freq.update(set(ids))
pres=[(g,freq[tokd[g]]) for g in drv["ensembl"] if g in tokd and tokd[g] in freq]
pres.sort(key=lambda x:-x[1]); cand=[g for g,_ in pres[:20]]
w(f"genes to perturb (expressed, top20 by cell-presence): {len(cand)}")
w("  e.g. "+", ".join(drv.loc[drv.ensembl.isin(cand[:8]),'symbol'].tolist()))

csm={"state_key":"aggr_state","start_state":"high","goal_state":"low","alt_states":[]}

# 3) state 嵌入字典（复用已存 pkl）
se_path=os.path.join(OUT,"state_embs.pkl")
if os.path.exists(se_path):
    state_embs=pickle.load(open(se_path,"rb")); w("reuse existing state_embs.pkl")
else:
    emb=EmbExtractor(model_type="Pretrained",num_classes=0,emb_mode="cls",
                     cell_emb_style="mean_pool",filter_data=None,max_ncells=None,
                     emb_layer=-1,forward_batch_size=10,nproc=2,
                     summary_stat="exact_mean", emb_label=None)
    w("computing state_embs_dict ...")
    state_embs=emb.get_state_embs(csm, MODEL, ds_path, OUT, "state_embs")

# 4) 锚定 perturbation：逐基因敲除（每次 filter 只需含该基因），仅作用于 high 态细胞子集
w("running per-gene in-silico deletion (high-state -> goal low) ...")
for i,g in enumerate(cand):
    try:
        isp=InSilicoPerturber(perturb_type="delete", genes_to_perturb=[g], combos=0,
                              model_type="Pretrained", num_classes=0, emb_mode="cls",
                              cell_emb_style="mean_pool",
                              cell_states_to_model=csm, state_embs_dict=state_embs,
                              max_ncells=200, emb_layer=-1, forward_batch_size=10, nproc=1)
        isp.perturb_data(MODEL, ds_path, OUT, f"anchor_{g}")
        w(f"  [{i+1}/{len(cand)}] {g} done")
    except Exception as e:
        w(f"  [{i+1}/{len(cand)}] {g} FAILED: {str(e)[:100]}")

# 5) 统计 goal-state shift 排序
st=InSilicoPerturberStats(mode="goal_state_shift", genes_perturbed=cand,
                          combos=0, cell_states_to_model=csm)
st.get_stats(OUT, None, OUT, "anchor_stats")
w("stats written -> "+OUT)

# 6) 对标经典 RF
try:
    sf=[f for f in os.listdir(OUT) if f.startswith("anchor_stats") and f.endswith(".csv")][0]
    s=pd.read_csv(os.path.join(OUT,sf))
    w("\n=== anchored shift stats head ===\n"+s.head(15).to_string())
    s.to_csv("results/deg/P2_anchored_shift_stats.csv",index=False)
    # merge with classical RF importance by ensembl/gene
    key=[c for c in s.columns if "Gene" in c or "gene" in c]
    w("\nstats columns: "+str(list(s.columns)))
    from scipy.stats import spearmanr
    m=s.merge(drv, left_on=key[0], right_on="ensembl", how="inner") if key else None
    if m is not None and len(m)>5:
        shiftcol=[c for c in s.columns if "Shift" in c or "shift" in c or "test" in c.lower()]
        if shiftcol:
            rho=spearmanr(m[shiftcol[0]], m["RF_importance"]).correlation
            w(f"\nFM goal-shift vs classical RF importance Spearman rho={rho:.3f} (concordance of FM & classical drivers)")
except Exception as e:
    w("benchmark step note: "+str(e)[:200])

w("\nALL-DONE")
log.close()
