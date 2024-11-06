import streamlit as st

from pathlib import Path
import streamlit as st

def read_markdown_file(markdown_file):
    return Path(markdown_file).read_text()

with open('bright-stars-select.tar.gz', 'rb') as f:
    st.download_button('Download this notebook and files',f)

intro_markdown = read_markdown_file("Bright-Stars-Select.md")
st.markdown(intro_markdown, unsafe_allow_html=True)
