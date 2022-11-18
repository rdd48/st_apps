import math
import subprocess
import pickle
import streamlit as st

st.set_page_config(
    page_title='Local EFI',
    page_icon=':atom_symbol:',
    layout='wide'
)

global blast_path 
blast_path = '/Users/robbydivine/opt/miniconda3/envs/lefi/bin/'

def process_fasta(fasta_file):

    fasta_dict = {}

    lines = fasta_file.getvalue().decode('utf-8').split('\n')

    for i in range(len(lines)):
        if lines[i][0] == '>':
            name = lines[i][1:]
            seq = ''
            new_index = i + 1
            while lines[new_index][0] != '>':
                seq += lines[new_index]
                new_index += 1
                if new_index == len(lines):
                    break
            
            name_only = name.split()[0]
            if '[' in name:
                species = name.split('[')[-1].replace(']', '').strip()
            else:
                species = ''

            if len(name.split()) > 1:
                descript = ' '.join(name.split()[1:])
                descript = descript.replace(f'[{species}]', '').strip()
            else:
                descript = ''

            fasta_dict[name_only] = (seq.replace('\n', ''), species, descript)

    return fasta_dict

def fasta_to_dicts(network_name):

    # keep info for xgmml
    nodes, edges = {}, {}

    # subprocess.run(f'{blast_path}makeblastdb -in {network_name} -dbtype prot'.split())
    out = subprocess.check_output(f"{blast_path}blastp -query {network_name} -subject {network_name} -matrix BLOSUM62 -outfmt '6 qseqid qlen sseqid slen bitscore length pident'", shell=True)

    for line in out.decode('utf-8').split('\n'):
        if line:
            n1, l1, n2, l2, bitscore, align_len, pident = line.split('\t')

            # nodes dict structure = name: (length,  sequence, species, descript)
            if n1 not in nodes:
                seq, species, descript = fasta_dict[n1]
                nodes[n1] = (l1, seq, species, descript)
            if n2 not in nodes:
                seq, species, descript = fasta_dict[n2]
                nodes[n2] = (l2, seq, species, descript)

            if n1 != n2:

                # https://github.com/EnzymeFunctionInitiative/EFITools/blob/master/lib/EFI/Util/AlignmentScore.pm
                # is formula from ^ actually?:
                align_score = -(math.log(l1 * l2) / math.log(10)) + (bitscore * math.log(2) / math.log(10)) # log must be ln

                # in perl:
                #    -(log($qlen * $slen) / log(10))
                #         +
                #     $bitscore * log(2) / log(10)
                # );

                # check for case where bitscore is so high that 2 ** -bitscore evaluates to 0.0
                # if (2 ** (-float(bitscore))):
                #     align_score = -math.log10((2. ** (-float(bitscore))) * float(l1) * float(l2))
                # else:
                #     align_score = -math.log10((2. ** (-float(1000.))) * float(l1) * float(l2))
                
                
                # edges dict structure = (name 1, name 2): (percent id, align_score, align_len)
                if (n1, n2) and (n2, n1) not in edges:
                    edges[(n1, n2)] = (pident, align_score, align_len)
                

    
    return nodes, edges

def dataset_analysis(nodes, edges):

    # keep info for descriptive graphs
    num_seqs, score_by_len, score_by_id, score_by_edge = {}, {}, {}, {}

    for v in nodes.values():
        seqlen = int(v[0])
        if seqlen not in num_seqs:
            num_seqs[seqlen] = 1
        else:
            num_seqs[seqlen] += 1
    
    for v in edges.values():
        pident, align_score, align_len = v

        pident = float(pident)
        align_score = int(align_score)
        align_len = int(align_len)

        align_score = round(align_score)

        if align_score not in score_by_len:
            score_by_len[align_score] = [align_len]
        else:
            v_copy = score_by_len[align_score].copy()
            v_copy.append(align_len)
            score_by_len[align_score] = v_copy

        if align_score not in score_by_id:
            score_by_id[align_score] = [pident]
        else:
            v_copy = score_by_id[align_score].copy()
            v_copy.append(pident)
            score_by_id[align_score] = v_copy

        if align_score not in score_by_edge:
            score_by_edge[align_score] = 1
        else:
            score_by_edge[align_score] += 1
    
    return num_seqs, score_by_len, score_by_id, score_by_edge

def pickle_dict(d, fname):
    with open(fname, 'wb') as fout:
        pickle.dump(d, fout)

def write_all_dicts(network_name, nodes, edges, num_seqs, score_by_len, score_by_id, score_by_edge):

    # # this would also work, but would create multiple files. 
    # pickle_dict(nodes, f'{network_name}_nodes.pickle')
    # pickle_dict(edges, f'{network_name}_edges.pickle')

    dict_of_dicts = {
        'nodes': nodes,
        'edges': edges,
        'num_seqs': num_seqs,
        'score_by_len': score_by_len,
        'score_by_id': score_by_id,
        'score_by_edge': score_by_edge
    }

    pickle_dict(dict_of_dicts, f'{network_name}_data.pickle')

def write_xgmml(network_name, nodes, edges, score_cutoff):

    with open(f'out/{network_name}_network.xgmml', 'w') as fout:

        # write header
        fout.write(f'<!-- Database: 1 -->\n\n<graph label="{network_name} Full Network" xmlns="http://www.cs.rpi.edu/XGMML">\n')

        # write nodes
        for node, v in nodes.items():
            length, seq, species, desc = v
            fout.write(f'  <node id="{node}" label="{node}">\n')
            # fout.write('    <att name="Sequence Source" type="string" value="USER" />\n')
            fout.write(f'    <att name="Sequence Length" type="integer" value="{length}" />\n')
            fout.write(f'    <att name="Species" type="string" value="{species}" />\n')
            fout.write(f'    <att name="Sequence" type="string" value="{seq}" />\n')
            fout.write('    <att type="list" name="Description">\n')
            fout.write(f'      <att type="string" name="Description" value="{desc}" />\n')
            fout.write('    </att>\n')
            fout.write('  </node>\n')
        
        # write edges
        for edge, v in edges.items():
            e1, e2 = edge
            pident, align_score, align_len = v
            if align_score > score_cutoff:
                fout.write(f'  <edge source="{e1}" target="{e2}" label="{e1},{e2}" id="{e1},{e2}">\n')
                fout.write(f'    <att name="%id" type="real" value="{pident}" />\n')
                fout.write(f'    <att name="alignment_score" type="real" value="{align_score}" />\n')
                fout.write(f'    <att name="alignment_len" type="real" value="{align_len}" />\n')
                fout.write('  </edge>\n')
        
        fout.write('</graph>')

st.markdown('WIP: generate .pickle files on streamlit.')
st.markdown('Current status: use https://github.com/rdd48/st_apps/blob/main/lefi/efi_cacao.py but update the blast_path on line 9')
st.markdown('Usage of the efi_cacao.py script: python efi_cacao.py PATH_TO_FASTA_INPUT')

# fasta_file = st.file_uploader('Upload the .fasta output from the initial dataset generation.')

# if fasta_file is not None:

#     fasta_dict = process_fasta(fasta_file)
    
#     network_name_full = fasta_file.name.strip()
#     network_name_short = fasta_file.name.split('/')[-1].split('.fasta')[0].strip()

#     nodes, edges = fasta_to_dicts(network_name_full)
#     num_seqs, score_by_len, score_by_id, score_by_edge = dataset_analysis(nodes, edges)

    # write_xgmml(network_name, nodes, edges, score_cutoff=90)
    # write_all_dicts(network_name, nodes, edges, num_seqs, score_by_len, score_by_id, score_by_edge)