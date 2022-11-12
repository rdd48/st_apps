import sys
import math
import subprocess

fasta_file = sys.argv[1]

global blast_path 
blast_path = '/Users/rddivine/opt/miniconda3/envs/bp/bin/'

def process_fasta(fasta_file):
    with open(fasta_file) as f:

        fasta_dict = {}

        lines = f.readlines()
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

def fasta_to_dicts(input_fasta):

    # keep info for xgmml
    nodes, edges = {}, {}

    fasta_dict = process_fasta(input_fasta)

    subprocess.run(f'{blast_path}makeblastdb -in {input_fasta} -dbtype prot'.split())
    out = subprocess.check_output(f"{blast_path}blastp -query {input_fasta} -db {input_fasta} -matrix BLOSUM62 -outfmt '6 qseqid qlen sseqid slen bitscore length pident'", shell=True)

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

                # check for case where bitscore is so high that 2 ** -bitscore evaluates to 0.0
                if (2 ** (-float(bitscore))):
                    align_score = -math.log10((2. ** (-float(bitscore))) * float(l1) * float(l2))
                else:
                    align_score = -math.log10((2. ** (-float(1000.))) * float(l1) * float(l2))
                
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


def write_dicts(network_name, nodes, edges, num_seqs, score_by_len, score_by_id, score_by_edge):

    with open(f'out/{network_name}_nodes.py', 'w') as fout:
        fout.write('nodes = {\n')
        for k,v in nodes.items():
            fout.write(f'    {k}: {v},\n')
        fout.write('}\n\n')

    with open(f'out/{network_name}_edges.py', 'w') as fout:
        fout.write('edges = {\n')
        for k,v in edges.items():
            fout.write(f'    {k}: {v},\n')
        fout.write('}\n\n')

    with open(f'out/{network_name}_data.py', 'w') as fout:
        fout.write('num_seqs = {\n')
        for k,v in num_seqs.items():
            fout.write(f'    {k}: {v},\n')
        fout.write('}\n\n')

        fout.write('score_by_len = {\n')
        for k,v in score_by_len.items():
            fout.write(f'    {k}: {v},\n')
        fout.write('}\n\n')

        fout.write('score_by_id = {\n')
        for k,v in score_by_id.items():
            fout.write(f'    {k}: {v},\n')
        fout.write('}\n\n')

        fout.write('score_by_edge = {\n')
        for k,v in score_by_edge.items():
            fout.write(f'    {k}: {v},\n')
        fout.write('}\n\n')



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


if __name__ == '__main__':
    nodes, edges = fasta_to_dicts(fasta_file)
    network_name = fasta_file.split('.fasta')[0].strip()
    num_seqs, score_by_len, score_by_id, score_by_edge = dataset_analysis(nodes, edges)

    write_xgmml(network_name, nodes, edges, score_cutoff=75)
    write_dicts(network_name, nodes, edges, num_seqs, score_by_len, score_by_id, score_by_edge)