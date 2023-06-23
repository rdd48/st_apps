from pytube import YouTube
import os


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


files = ['https://www.youtube.com/watch?v=CB2WJG-YvP4']
for f in files:
    download_360p_mp4_videos(f)

