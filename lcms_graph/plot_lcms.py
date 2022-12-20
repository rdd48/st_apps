import streamlit as st
import pandas as pd
import plotly.express as px

# try:
#     import plotly.express as px
# except ImportError:
#     subprocess.check_call([sys.executable, "-m", "pip", "install", "plotly"])
#     subprocess.check_call([sys.executable, "-m", "pip", "install", "xlrd"])
# finally:
#     import plotly.express as px


st.set_page_config(
    page_title='LC/MS Plotter',
    page_icon=':chart_with_upwards_trend:',
    layout='wide'
)

# session states
if 'y_options' not in st.session_state:
    st.session_state.y_options = None

def set_all_y_options(all_options):
    st.session_state.y_options = all_options

def set_no_y_options():
    st.session_state.y_options = []

# process input to df
@st.cache
def csv_to_df(csv_file):
    df = pd.read_csv(
        csv_file,
        names = ['x', csv_file.name.replace('.csv', '')]
    )

    # drop extra l columns
    # df.drop(list(df.filter(regex = 'l\.')), axis = 1, inplace = True)
    return df

st.title(':chart_with_upwards_trend: LC/MS Plotter')

uploaded_csvs = st.file_uploader(
    'Upload your csv file here:',
    type='csv',
    accept_multiple_files=True
)


if uploaded_csvs:
    dfs = []
    for csv in uploaded_csvs:
        dfs.append(csv_to_df(csv))
    
    df = dfs[0]
    if len(dfs) > 1:
        for df_to_merge in dfs:
            df = pd.merge(df.copy(), df_to_merge, 'outer')

    st.sidebar.header('Please filter here:')
    y_options = st.sidebar.multiselect(
        'Select the data to plot:',
        options = df.columns.to_list()[1:],
        default = df.columns.to_list()[1]
    )

    line_chart = px.line(
        df,
        x='x',
        y=y_options,
    )

    line_chart.update_layout(width=1000)
    #line_chart.update_layout(autosize=True)

    # all_button = st.sidebar.button('Select all', on_click=set_all_y_options, args=(df.columns.to_list()[1:]))
    # remove_button = st.sidebar.button('Remove all', on_click=set_no_y_options)

    st.plotly_chart(line_chart)