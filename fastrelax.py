import streamlit as st
import py3Dmol
from stmol import showmol

import pyrosetta
from pyrosetta.rosetta.protocols.relax import FastRelax

def pose_to_str(pose):
    ss = pyrosetta.rosetta.std.stringstream()
    pyrosetta.rosetta.core.io.pdb.dump_pdb( pose, ss )
    return ss.str()

def display_pdb(pdb_str):
    view = py3Dmol.view(width=400, height=400)
    view.addModelsAsFrames(pdb_str)
    view.setStyle({'cartoon':{'color':'spectrum'}})
    view.zoomTo()
    view.show()

    showmol(view, width=400, height=400)

def run_fastrelax():
    if pdb_code == '':
        st.markdown('Please provide an RCSB pdb code above!')
        return

    pyrosetta.init()

    sfxn = pyrosetta.get_fa_scorefxn()
    fr = FastRelax()
    fr.set_scorefxn(sfxn)

    # pose = pyrosetta.pose_from_pdb('1ubq.pdb')
    ## TO DO: write my own fxn here to get the data from RCSB
    ## see https://github.com/RosettaCommons/main/blob/master/source/src/python/PyRosetta/src/pyrosetta/toolbox/rcsb.py
    ## especially the download_from_web and pose_from_rcsb fxns
    ## pdb_url = "http://files.rcsb.org/download/" + pdb_code + ".pdb"
    pose_from_rcsb = pyrosetta.toolbox.rcsb.pose_from_rcsb(pdb_code)

    pdb_str = pose_to_str(pose_from_rcsb)
    display_pdb(pdb_str)
    
    fr.apply(pose_from_rcsb)   

    pdb_str_relax = pose_to_str(pose_from_rcsb)
    display_pdb(pdb_str_relax)


pdb_code = st.text_input('Enter 4-character PDB code:', value='')
run_fr_button = st.button('Click to FastRelax', on_click=run_fastrelax)


