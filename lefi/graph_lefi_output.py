import streamlit as st
import pickle
import plotly.express as px
import plotly.graph_objects as go

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

uploaded_pickle = st.file_uploader('Upload the .pickle output from the initial dataset generation.')

if uploaded_pickle is not None:
    data = pickle.load(uploaded_pickle)

    # bin_size = st.sidebar.slider(
    #     'Select the number of bins for the histograms', 
    #     min_value=10, 
    #     max_value=max(data['num_seqs'].keys()) // 2,
    #     value=50 
    # )

    bin_size = st.sidebar.number_input(
        'Select the number of bins for the histograms', 
        min_value=10, 
        max_value=max(data['num_seqs'].keys()) // 2,
        value=50 
    )

    resized_ns = resize_by_num_boxes(data['num_seqs'], bin_size)
    seq_len_hist = px.bar(
        x = resized_ns.keys(),
        y = resized_ns.values(),
        # x = data['num_seqs'].keys(),
        # y = data['num_seqs'].values(),
        labels={'x':'Sequence Length', 'y':'Count'})

    score_len_hist = go.Figure()
    resized_sbl = resize_by_num_boxes_boxplot(data['score_by_len'], bin_size)
    for k in sorted(resized_sbl.keys()):
        score_len_hist.add_trace(go.Box(
            y=resized_sbl[k],
            name=k,
            marker_color='royalblue',
        ))
    # for k in sorted(data['score_by_len'].keys()):
    #     score_len_hist.add_trace(go.Box(
    #         y=data['score_by_len'][k],
    #         name=k,
    #         marker_color='royalblue',
    #     ))
    score_len_hist.update_layout(
        xaxis={'title': 'Alignment Score'},
        yaxis={'title': 'Alignment Length'},
        showlegend=False
    )

    score_id_hist = go.Figure()
    resized_sbi = resize_by_num_boxes_boxplot(data['score_by_id'], bin_size)
    for k in sorted(resized_sbi.keys()):
        score_id_hist.add_trace(go.Box(
            y=resized_sbi[k],
            name=k,
            marker_color='royalblue',
        ))
    # for k in sorted(data['score_by_id'].keys()):
    #     score_id_hist.add_trace(go.Box(
    #         y=data['score_by_id'][k],
    #         name=k,
    #         marker_color='royalblue',
    #     ))
    score_id_hist.update_layout(
        xaxis={'title': 'Alignment Score'},
        yaxis={'title': 'Percent Identity'},
        showlegend=False
    )

    
    reverse_scores = sorted(data['score_by_edge'].keys(), reverse=True)
    y_edges = []
    curr_edge_num = 0
    for s in reverse_scores:
        curr_edge_num += data['score_by_edge'][s]
        y_edges.append(curr_edge_num)

    align_score_line = px.line(
        x = reverse_scores,
        y = y_edges,
        labels={'x':'Alignment Score', 'y':'Number of Edges'})

    resized_sbe = resize_by_num_boxes(data['score_by_edge'], bin_size)
    edge_count_hist = px.bar(
        x = resized_sbe.keys(),
        y = resized_sbe.values(),
        labels={'x':'Alignment Score', 'y':'Number of Edges'})


    st.plotly_chart(seq_len_hist)
    st.plotly_chart(score_len_hist)
    st.plotly_chart(score_id_hist)
    st.plotly_chart(align_score_line)
    st.plotly_chart(edge_count_hist)