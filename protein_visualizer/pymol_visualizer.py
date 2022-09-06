from pdb import Pdb
import streamlit as st
import py3Dmol
from stmol import showmol
from urllib import request
import IPython
# import gzip

# local files
from utils import utils
from pdb_info import pdb_names_dict

# TODO:
# * save/render images -- maybe through pymol api?
# * save/render gifs
# * clean up code
# * host 

st.set_page_config(
    page_title='Pymol Visualizer',
    page_icon=':atom_symbol:',
    layout='wide'
)

if 'pdb_state' not in st.session_state:
    st.session_state.pdb_state = None

def set_pdb_state(code):
    if code in ('rcsb_code', 'user_input', 'rcsb_desc'):
        st.session_state.pdb_state = code
    else:
        st.session_state.pdb_state = None

# st.markdown('#### :atom_symbol:  Pymol Visualizer  :atom_symbol:')
mobile_size = st.checkbox('Optimize size for mobile?', value=False)

# globals
width = 250 if mobile_size else 600
height = 250 if mobile_size else 400

col1, col2, col3 = st.columns([1,1,1])

# col1, col2 = st.columns([1,2])
with col1:
    defaults_to_pdb = {
        'Top7': '1QYS',
        'Kumamolisin': '1SIO',
        'Designed retroaldolase': '3HOJ'
    }
    selected_pdb = st.selectbox('Pick from common options:', options=['', *defaults_to_pdb.keys()], on_change=set_pdb_state, args=['user_input'])

with col2:
    rcsb_input = st.text_input('Or enter PDB code:', on_change=set_pdb_state, args=['rcsb_code'])
    rcsb_input = rcsb_input.strip()

with col3:
    rcsb_desc = st.selectbox('Or search by protein name (slightly laggy)', options=[*pdb_names_dict.all_pdbs_dict.keys()], on_change=set_pdb_state, args=['rcsb_desc'])

if st.session_state.pdb_state == 'rcsb_code':
    pdb = rcsb_input.strip().upper()
elif st.session_state.pdb_state == 'user_input' and selected_pdb:
    pdb = defaults_to_pdb[selected_pdb].strip().upper()
elif st.session_state.pdb_state == 'rcsb_desc':
    pdb = pdb_names_dict.all_pdbs_dict[rcsb_desc]
else:
    pdb = None

st.sidebar.header('Select options:')
color_selector = st.sidebar.selectbox('Select color from picker or from pymol defaults:', ('Spectrum', 'Color by Chain', 'Color Picker',))
spinning = st.sidebar.selectbox('Set model spin', ('Off', 'On'))

if pdb and pdb in pdb_names_dict.all_pdbs_dict.values():
    url = f'https://files.rcsb.org/download/{pdb}.pdb1.gz'
    r = request.urlopen(url=url)

    PdbObj = utils.parse_pdb(r, pdb)

    view = py3Dmol.view(width=width, height=height)
    view.addModelsAsFrames(PdbObj.display_pdb_lines())
    # view.setStyle({'model': -1}, {'cartoon': {'color': 'spectrum'}})

    # display pdb element title
    st.markdown(f'**[{pdb.upper()}](https://www.rcsb.org/structure/{pdb.upper()}): {utils.get_name(pdb)}**')

    # PdbObj = utils.parse_pdb(pdb, pdb_lines)
    chains = PdbObj.get_chains()

    def color_model(color):
        all_backbone = PdbObj.backbone_atoms
        for atom in all_backbone:
            view.setStyle({'model': -1, 'serial': atom}, {'cartoon': {'color': color}})
        all_hetatm = PdbObj.het_atoms.values()
        for idx, atomlist in enumerate(all_hetatm):
            for atom in atomlist:
                view.setStyle({'model': -1, 'serial': atom}, {'stick': {'color': sidechain_colors[idx]}})


    colors = ['green', 'cyan', 'magenta', 'yellow', 'wheat', 'purple', 'grey', 'lightpink', 'blue', 'red']
    sidechain_colors = ['red', 'cyan', 'purple', 'yellow', 'wheat', 'green', 'blue', 'lightpink', 'cyan', 'magenta', ]
    if color_selector == 'Color Picker':
        color = st.sidebar.color_picker('Select color from picker:', value='#00f900')
        color_model(color)
        # view.setStyle({'model': -1}, {'cartoon': {'color': color}})
    elif color_selector == 'Color by Chain':
        for idx, c in enumerate(chains):
            view.setStyle({'chain': c}, {'cartoon': {'color': colors[idx % len(colors)]}})
        all_hetatm = PdbObj.het_atoms.values()
        for idx, atomlist in enumerate(all_hetatm):
            for atom in atomlist:
                view.setStyle({'model': -1, 'serial': atom}, {'stick': {'color': sidechain_colors[idx]}})
    elif color_selector == 'Spectrum':
        color = color_selector.lower()
        color_model(color)
        # view.setStyle({'model': -1}, {'cartoon': {'color': color}})
    

    if spinning == 'On':
        spin_axis = st.sidebar.selectbox('Which axis?', ('x', 'y', 'z'))
        speed = st.sidebar.slider('How fast? (Default is 1)', min_value = 0.25, max_value=2., value=1., step=0.25)
        spin_reverse = st.sidebar.selectbox('Reverse?', ('False', 'True'))
        spin_direction = 1.0 if spin_reverse == 'False' else -1.0
        view.spin(spin_axis, speed * spin_direction)

    button = st.sidebar.button('Reset view', on_click=view.zoomTo)

    view.zoomTo()
    view.show()

    script = '''<iframe>
            <script>
            var png = viewer_{0}.pngURI()
            $('iframe').attr('src', png)
            </script>'''.format(view.uniqueid)

    IPython.display.publish_display_data({'application/3dmoljs_load.v0':script, 'text/html': script},metadata={})
    

    
    showmol(view, width=width, height=height)

    st.download_button('Download sequence as fasta', PdbObj.get_sequence(), file_name=f'{pdb}.fasta')

    aa_options = sorted([aa for aa in PdbObj.sidechain_atoms_by_aas.keys() if aa not in ('UNK', 'HOH')])
    aas_highlight = st.multiselect('Select amino acids to show', aa_options)
    
    if aas_highlight:
        sidechain_atoms = {}
        for aa in aas_highlight:
            sidechain_atoms[aa] = PdbObj.get_sidechain_atoms_by_aa(aa)
            backbone_atoms = PdbObj.backbone_atoms_no_o2
            
        # for idx in range(PdbObj.num_atoms):
            # view.setStyle({'model': -1, 'serial': idx+1}, {'stick': {'color': 'red'}})
        view2 = py3Dmol.view(width=width, height=height)
        view2.addModelsAsFrames(PdbObj.parsed_pdb)
        for atom in backbone_atoms:
            view2.setStyle({'model': -1, 'serial': atom}, {'stick': {'color': color}})
        for idx, atom in enumerate(sidechain_atoms.values()):
            # view.setStyle({'model': -1, 'serial': idx+1}, {'stick': {'color': 'black'}})
            view2.setStyle({'model': -1, 'serial': atom}, {'stick': {'color': sidechain_colors[idx]}})
        
        view2.zoomTo()
        view2.show()
        showmol(view2, width=width, height=height)



# center the title
# st.markdown('<style>.e16nr0p30 {text-align: center;}</style>', unsafe_allow_html=True)

# adjust the view size and add a border
st.markdown('<style>iframe {border: solid; display: block; margin-left: auto; margin-right: auto; padding-left: 0px;}</style>', unsafe_allow_html=True)