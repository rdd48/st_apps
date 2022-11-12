import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


from out import cdma_blast_cdhit_data as data

# num_seqs, score_by_len, score_by_id, score_by_edge

# df = pd.DataFrame(data.num_seqs)
# df.head()

seq_len_hist = px.bar(
    # df,
    x = data.num_seqs.keys(),
    y = data.num_seqs.values(),
    labels={'x':'Sequence Length', 'y':'Count'})

score_len_hist = go.Figure()
for k in sorted(data.score_by_len.keys()):
    score_len_hist.add_trace(go.Box(
        y=data.score_by_len[k],
        name=k,
    ))
score_len_hist.update_layout(
    xaxis={'title': 'Alignment Score'},
    yaxis={'title': 'Alignment Length'}
)

score_id_hist = go.Figure()
for k in sorted(data.score_by_id.keys()):
    score_id_hist.add_trace(go.Box(
        y=data.score_by_id[k],
        name=k,
    ))
score_id_hist.update_layout(
    xaxis={'title': 'Alignment Score'},
    yaxis={'title': 'Percent Identity'}
)

reverse_scores = sorted(data.score_by_edge.keys(), reverse=True)
y_edges = []
curr_edge_num = 0
for s in reverse_scores:
    curr_edge_num += data.score_by_edge[s]
    y_edges.append(curr_edge_num)

align_score_line = px.line(
    x = reverse_scores,
    y = y_edges,
    labels={'x':'Alignment Score', 'y':'Number of Edges'})

edge_count_hist = px.bar(
    x = data.score_by_edge.keys(),
    y = data.score_by_edge.values(),
    labels={'x':'Alignment Score', 'y':'Number of Edges'})


st.plotly_chart(seq_len_hist)
st.plotly_chart(score_len_hist)
st.plotly_chart(score_id_hist)
st.plotly_chart(align_score_line)
st.plotly_chart(edge_count_hist)