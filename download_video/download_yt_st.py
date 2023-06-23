from pytube import YouTube
import streamlit as st
import os

st.set_page_config(
    page_title='Download Video',
    layout='wide'
)

def download_360p_mp4_videos(url: str, outpath: str = ''):
    if not outpath:
        outpath = os.getcwd()
    yt = YouTube(url)
    ys = yt.streams.get_highest_resolution()
    ys.download()
    # yt = YouTube(url)
    # mp4files = yt.filter('mp4') 
    # yt.set_filename('GeeksforGeeks Video')
    # d_video = yt.get(mp4files[-1].extension,mp4files[-1].resolution)
    # d_video.download(outpath)  
    # yt.streams.filter(file_extension="mp4").get_by_resolution("360p").download(outpath)

def button_callback(files):
    files = [f.strip() for f in files.split('\n')]
    for f in files:
        print(f)
        download_360p_mp4_videos(f)

files = st.text_area('Enter videos to download here, each on its own line')
button = st.button(
    'Click to download videos',
    on_click=button_callback,
    args=[files]
)
