import streamlit as st
from pathlib import Path


st.set_page_config(
        page_title="AOdb",
        page_icon="images/Starcutout.png",
        layout="wide",
    )

sidebar_logo = 'images/Starcutout.png'
st.logo(sidebar_logo, size='large')

def read_markdown_file(markdown_file):
    return Path(markdown_file).read_text()

with open('bright-stars-select.tar.gz', 'rb') as f:
    st.download_button('Download this notebook and files',f)

intro_markdown = read_markdown_file("Bright-Stars-Select.md")
st.markdown(intro_markdown, unsafe_allow_html=True)
