#!/usr/bin/env python
# 17c — 下载 Geneformer 真·字典(gc104M) + V2-104M 权重，覆盖包内 LFS 指针
import os, shutil, geneformer
from huggingface_hub import hf_hub_download
REPO="ctheodoris/Geneformer"
pkg=os.path.dirname(geneformer.__file__)
print("pkg dir:",pkg)
dicts=["gene_median_dictionary_gc104M.pkl","token_dictionary_gc104M.pkl",
       "ensembl_mapping_dict_gc104M.pkl","gene_name_id_dict_gc104M.pkl"]
for d in dicts:
    p=hf_hub_download(REPO, filename=f"geneformer/{d}")
    shutil.copy(p, os.path.join(pkg,d))
    print("dict ok:",d, os.path.getsize(os.path.join(pkg,d)),"bytes")
# model V2-104M
mdir=os.path.join("results/scrna/gf_model_V2-104M"); os.makedirs(mdir,exist_ok=True)
for f in ["config.json","model.safetensors","training_args.bin"]:
    try:
        p=hf_hub_download(REPO, filename=f"Geneformer-V2-104M/{f}")
        shutil.copy(p, os.path.join(mdir,f)); print("model ok:",f, os.path.getsize(os.path.join(mdir,f)),"bytes")
    except Exception as e:
        print("model skip",f,str(e)[:80])
print("ASSETS-DONE", mdir)
