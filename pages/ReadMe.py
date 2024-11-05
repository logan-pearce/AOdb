import streamlit as st
import pandas as pd
import numpy as np

st.set_page_config(
        page_title="AOEngDB",
        page_icon="",
        layout="wide",
    )

st.title('Database of Bright Single Stars for AO Engineering')

st.markdown(
    """
    SLSdb is an interactive database of all known Sirius-Like Systems -- white dwarfs with non-interacting main sequence star companions of spectral type earlier than M. View the Tutorials page for how to interact with the db; the Derivation page walks through how we adapted multiple catalogs into the db.

### Contibuting
To contribute to slsdb please email Logan Pearce at lapearce@umich.edu

### Future Upgrades
 - Allow users to submit new SLS

"""
)

