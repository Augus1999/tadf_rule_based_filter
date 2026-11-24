# -*- coding: utf-8 -*-
# Author: Nianze TAO
"""
Stricter filtering out non-novel and unstable structures.
"""
import time
import requests
from contextlib import redirect_stderr
from io import StringIO
from concurrent.futures import ThreadPoolExecutor
from rdkit import rdBase
from rdkit.Chem import MolFromSmiles, AddHs, AllChem


rdBase.LogToPythonStderr()

sample_file = "YOUR/SAMPLE/FILE_filtered.csv"


def _embed(m, num_conformers=1):
    m = AddHs(m)
    with StringIO() as buf:
        with redirect_stderr(buf):
            AllChem.EmbedMultipleConfs(m, numConfs=num_conformers, randomSeed=42)
            for conf_id in range(num_conformers):
                ff = AllChem.UFFGetMoleculeForceField(m, confId=conf_id)
                ff.CalcEnergy()
            res = buf.getvalue()
            if "UFFTYPER" in res:
                return False
    return True


data_ = []
with open(sample_file, "r", encoding="utf-8") as f:
    data = f.readlines()
for i in data:
    smi = i.replace("\n", "").split(",")[0]
    if (
        a := requests.post(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/cids/txt",
            data={"smiles": smi},
        ).text.strip()
    ) != "0":
        print(smi, a)
        time.sleep(1.0)
        continue
    mol = MolFromSmiles(smi)
    with ThreadPoolExecutor() as excutor:
        future = excutor.submit(_embed, mol)
        try:
            r = future.result(timeout=30)
            if r:
                data_.append(i)
        except Exception:
            ...
with open(sample_file.replace("filtered", "post_filtered"), "w", encoding="utf-8") as g:
    g.write("\n".join(data_))
