import streamlit as st
from pathlib import Path


st.set_page_config(
        page_title="AOdb",
        page_icon="images/Starcutout.png",
        layout="wide",
    )

sidebar_logo = 'images/Starcutout.png'
st.logo(sidebar_logo, size='large')

''' # In Memoriam
### To those bright star catalog members we lost along the way.
'''

st.text("")
st.text("")
st.text("")
st.text("")
cols = st.columns(2)
with cols[0]:
    '''### E Hya: triple'''
    st.image('images/EHya.png', width = 650)

with cols[1]:
    '''### alf Eri: binary'''
    st.image('images/alfEri.png', width = 600)