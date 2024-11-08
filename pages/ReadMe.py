import streamlit as st
import pandas as pd
import numpy as np

st.set_page_config(
        page_title="AOdb",
        page_icon="images/Starcutout.png",
        layout="wide",
    )

sidebar_logo = 'images/Starcutout.png'
st.logo(sidebar_logo, size='large')

st.title('Database of Bright Single Stars for AO Engineering at Las Campanas')

st.markdown(
"""
This is an interactive compilation of bright stars spanning all RAs and DECs near overhead at LCO. The Derivation page explains how the stars were selected.  
The 'Num' column corresponds to the index that will be in the TCS catalog.  The plot below the db shows the stars, hover over each point to see information including index, name, RA, and magnitudes. Enter the desired filename in the text box then click "Generate TCS Catalog" to download the full catalog in the correct TCS format. ::NOTE:: Filename must end in '.cat' to be used by telescope TCS.

The "CreateTCSCat" page lets you create TCS catalogs of objects providing just the Simbad names.  Supply a list of Simbad-resolvable names and download a TCS catalog for those objects.

### Contibuting
To update the catalog please email Logan Pearce at lapearce@umich.edu

### Future Upgrades
 - Allow users to update the database in real time


"""
)

'''## Example basic SQL queries'''

''' You can use SQL to select only desired stars from the full catalog and generate a new TCS catalog from the subset.'''

''' #### Ex: Select only stars within a specific RA range '''
strg = "SELECT * FROM aodb WHERE ra < 180 and ra > 50"
st.code(strg, language="sql")
'''Note: in order to generate a TCS catalog from this sql query, you must select all columns (*)'''

"""Written by Logan Pearce, 2024"""