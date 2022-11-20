import streamlit as st
import pickle
# import plotly.express as px
# import plotly.graph_objects as go

st.set_page_config(
    page_title='Local EFI',
    page_icon=':atom_symbol:',
    layout='wide'
)

def resize_by_num_boxes(d, box_nums):
    max_xaxis = max(d.keys())
    if box_nums >= max_xaxis:
        return d

    box_size = max_xaxis // box_nums
    if box_size == 1:
        return d

    resized_d = {(k // box_size) * box_size: 0 for k in range(0, max_xaxis, box_size)}

    for k, v in d.items():
        resized_d[(k // box_size) * box_size] += v
    
    return resized_d

def resize_by_num_boxes_boxplot(d, box_nums):
    max_xaxis = max(d.keys())
    if box_nums >= max_xaxis:
        return d

    box_size = max_xaxis // box_nums
    if box_size == 1:
        return d

    resized_d = {(k // box_size) * box_size: [] for k in range(0, max_xaxis, box_size)}

    for k, v in d.items():
        new_k = (k // box_size) * box_size
        v_copy = resized_d[new_k].copy()
        v_copy += v
        resized_d[new_k] = v_copy
    
    return resized_d

@st.cache
def write_xgmml(network_name, nodes, edges, score_cutoff):

    # write header
    fout = f'<!-- Database: 1 -->\n\n<graph label="{network_name} Full Network" xmlns="http://www.cs.rpi.edu/XGMML">\n'

    # write nodes
    for node, v in nodes.items():
        length, seq, species, desc = v
        fout += f'  <node id="{node}" label="{node}">\n'
        # fout += '    <att name="Sequence Source" type="string" value="USER" />\n'
        fout += f'    <att name="Sequence Length" type="integer" value="{length}" />\n'
        fout += f'    <att name="Species" type="string" value="{species}" />\n'
        fout += f'    <att name="Sequence" type="string" value="{seq}" />\n'
        fout += '    <att type="list" name="Description">\n'
        fout += f'      <att type="string" name="Description" value="{desc}" />\n'
        fout += '    </att>\n'
        fout += '  </node>\n'
    
    # write edges
    num_edges = 0
    for edge, v in edges.items():
        e1, e2 = edge
        pident, align_score, align_len = v
        if align_score > score_cutoff:
            num_edges += 1
            fout += f'  <edge source="{e1}" target="{e2}" label="{e1},{e2}" id="{e1},{e2}">\n'
            fout += f'    <att name="%id" type="real" value="{pident}" />\n'
            fout += f'    <att name="alignment_score" type="real" value="{align_score}" />\n'
            fout += f'    <att name="alignment_len" type="real" value="{align_len}" />\n'
            fout += '  </edge>\n'
    
    fout += '</graph>'

    return fout, num_edges

uploaded_pickle = st.file_uploader('Upload the .pickle output from the initial dataset generation.')

# if uploaded_pickle is not None:
#     data = pickle.load(uploaded_pickle)

#     sorted_scores = sorted(data['score_by_len'].keys())
#     if len(sorted_scores) > 1:
#         score_cutoff = st.slider(
#             'Set the score cutoff for each edge', 
#             min_value=sorted_scores[0], 
#             max_value=sorted_scores[-1],
#             value=sorted_scores[len(sorted_scores) // 2]
#             )
#     else:
#         score_cutoff = sorted_scores[0]
#         st.markdown(f'The only score is {score_cutoff}!')

#     network_name = uploaded_pickle.name.replace('.pickle', '.xgmml')
#     fout, num_edges = write_xgmml(network_name, data['nodes'], data['edges'], score_cutoff)

#     st.markdown(f'**{num_edges}** edges using a score cutoff of {score_cutoff}')
#     st.download_button(
#         'Download .xgmml for Cytoscape graphing',
#         fout,
#         file_name=network_name
#         )

#     max_bins = max(data['num_seqs'].keys())
#     bin_size = st.sidebar.number_input(
#         'Select the number of bins for the histograms', 
#         min_value=10 if max_bins // 5 > 10 else 1, 
#         max_value=max_bins // 2,
#         value=max_bins // 5
#     )

#     resized_ns = resize_by_num_boxes(data['num_seqs'], bin_size)
#     seq_len_hist = px.bar(
#         x = resized_ns.keys(),
#         y = resized_ns.values(),
#         labels={'x':'Sequence Length', 'y':'Count'})

#     score_len_hist = go.Figure()
#     resized_sbl = resize_by_num_boxes_boxplot(data['score_by_len'], bin_size)
#     for k in sorted(resized_sbl.keys()):
#         marker_color = 'lightblue' if k < score_cutoff else 'royalblue'
#         score_len_hist.add_trace(go.Box(
#             y=resized_sbl[k],
#             name=k,
#             marker_color=marker_color,
#         ))
#     score_len_hist.update_layout(
#         xaxis={'title': 'Alignment Score'},
#         yaxis={'title': 'Alignment Length'},
#         showlegend=False
#     )

#     score_id_hist = go.Figure()
#     resized_sbi = resize_by_num_boxes_boxplot(data['score_by_id'], bin_size)
#     for k in sorted(resized_sbi.keys()):
#         marker_color = 'lightblue' if k < score_cutoff else 'royalblue'
#         score_id_hist.add_trace(go.Box(
#             y=resized_sbi[k],
#             name=k,
#             marker_color=marker_color,
#         ))
#     score_id_hist.update_layout(
#         xaxis={'title': 'Alignment Score'},
#         yaxis={'title': 'Percent Identity'},
#         showlegend=False
#     )

#     reverse_scores = sorted(data['score_by_edge'].keys(), reverse=True)
#     y_edges = []
#     curr_edge_num = 0
#     for s in reverse_scores:
#         curr_edge_num += data['score_by_edge'][s]
#         y_edges.append(curr_edge_num)

#     align_score_line = px.line(
#         x = reverse_scores,
#         y = y_edges,
#         labels={'x':'Alignment Score', 'y':'Number of Edges (Cumulative)'},
#         # marker_color=num_edge_colors
#     )

#     edge_count_hist = go.Figure()
#     resized_sbe = resize_by_num_boxes(data['score_by_edge'], bin_size)
#     for k in sorted(resized_sbe.keys()):
#         marker_color = 'lightblue' if k < score_cutoff else 'royalblue'
#         edge_count_hist.add_trace(go.Bar(
#             x=[k], # these have to be lists? dumb
#             y=[resized_sbe[k]], # these have to be lists? dumb
#             name=k,
#             marker_color=marker_color,
#         ))
#     edge_count_hist.update_layout(
#         xaxis={'title': 'Alignment Score'},
#         yaxis={'title': 'Number of Edges (per bin)'},
#         showlegend=False
#     )

    # st.markdown('Node dataset analysis')
    # st.plotly_chart(seq_len_hist)

    # st.markdown('Edge dataset analysis')
    # st.plotly_chart(score_len_hist)
    # st.plotly_chart(score_id_hist)
    
    # display_vline = st.checkbox('Display score cutoff on below graph?', value=True)
    # if display_vline:
    #     align_score_line.add_vline(x=score_cutoff)
    # st.plotly_chart(align_score_line)

    # st.plotly_chart(edge_count_hist)


