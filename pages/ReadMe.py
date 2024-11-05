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
    This is an interactive compilation of bright stars spanning all RAs and DECs near overhead at LCO. The Derivation page explains how the stars were selected.  
    The 'Num' column corresponds to the index that will be in the TCS catalog.  The plot below the db shows the stars, hover over each point to see information including index, name, RA, and magnitudes.  Click "Generate TCS Catalog" to download the full catalog in the correct TCS format.

    The "CreateTCSCat" page lets you create TCS catalogs of objects providing just the Simbad names.  Supply a list of Simbad-resolvable names and download a TCS catalog for those objects.

### Contibuting
To update the catalog please email Logan Pearce at lapearce@umich.edu

### Future Upgrades
 - Allow SQL queries of database to downselect and produce a new TCS catalog
 - Allow users to update the database in real time


"""
)

