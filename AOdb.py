import streamlit as st
import pandas as pd
import numpy as np
from bokeh.plotting import figure
from bokeh.io import curdoc
from bokeh.palettes import Magma256, Viridis256
import sqlite3
from streamlit import session_state
import astropy.units as u

st.set_page_config(
        page_title="AOdb",
        page_icon="",
        layout="wide",
    )

st.title('Database of Bright Single Stars for AO Engineering from Las Campanas')

import sqlite3
import pandas as pd

aodb = pd.read_csv('Bright-AO-Stars.csv')

#aodb = pd.read_csv('aodb.csv')


import sqlite3
from sqlalchemy import create_engine, text
conn = sqlite3.connect('aodb.db')
engine = create_engine("sqlite:///aodb.db")
aodb.to_sql(name = 'aodb', con=engine, index=False, if_exists='replace')


def querySQL(string):
    with engine.connect() as conn:
        result = conn.execute(text(string)).fetchall()
        session_state['db'] = pd.DataFrame(result)
        st.dataframe(session_state['db'])

st.text_input(r"$\textsf{\Large SQL Query String}$", key='sqlquerystring')

session_state['db'] = aodb

#session_state
if session_state['sqlquerystring'] == '':
    session_state['db'] = aodb
    st.dataframe(session_state['db'])
else:
    with st.form(key="aodbsql"):
        st.form_submit_button('Query', on_click=querySQL(session_state['sqlquerystring']))




def GenerateCat(dataframe):

    pdcat = dataframe.copy()
    #pdcat['RA deg'],pdcat['Dec deg'] = np.nan,np.nan
    pdcat['pmra s/yr'], pdcat['pmdec arcsec/yr'] = np.nan, np.nan

    from astropy.coordinates import Angle
    ## convert pmra in mas/yr into s/yr and pmdec in mas/yr to arcsec/yr:
    # For each object:
    for i in range(len(pdcat)):
        # Create an astropy angl object:
        a = Angle(pdcat['pmra'][i],u.mas)
        # Convert to hms:
        a2 = a.hms
        # add up the seconds (a2[0] and a2[1] are most likely 0 but just in case):
        a3 = a2[0]*u.hr.to(u.s) + a2[1]*u.min.to(u.s) + a2[2]
        # put into table:
        pdcat.loc[i,'pmra s/yr'] = a3
        
        # Dec is easier:
        a = pdcat.loc[i,'pmdec']*u.mas.to(u.arcsec)
        # put into table:
        pdcat.loc[i,'pmdec arcsec/yr'] = a

    pdcat['num'] = np.arange(1,len(pdcat)+1,1)
    pdcat['Name'] = pdcat['simbad_name']

    pdcat_out = pdcat[['num','Name']]
    pdcat_out['ra'] = pdcat['ra_hms']
    pdcat_out['dec'] = pdcat['dec_dms']
    pdcat_out['Equinox'] = '2000'
    pdcat_out['pmra'] = pdcat['pmra s/yr']
    pdcat_out['pmdec'] = pdcat['pmdec arcsec/yr'] 
    pdcat_out['rotang'] = '0'
    pdcat_out['rot_mode'] = 'GRV'
    pdcat_out['RA_probe1'],pdcat_out['Dec_probe1'] = '00:00:00.00',  '+00:00:00.0'
    pdcat_out['equinox'] = '2000'
    pdcat_out['RA_probe2'],pdcat_out['Dec_probe2'] = '00:00:00.00',  '+00:00:00.0'
    pdcat_out['equinox '] = '2000'
    pdcat_out['epoch'] = '2000'
    pdcat_out = pdcat_out.to_csv(index=False, sep='\t')
    return pdcat_out


filename = st.text_input('TCS Cat Filename', 'Bright-AO-Stars.cat')
st.download_button(label='Generate TCS catalog', data=GenerateCat(session_state['db']), file_name=filename)




##################################### RA/DEC Plot::::
from bokeh.models import LinearColorMapper, ColumnDataSource, LinearInterpolator, ColorBar, Label
from bokeh.transform import linear_cmap, log_cmap
multiplier = 100
datadf = pd.DataFrame(data={'plotx':session_state['db']['ra'], 
                        'ploty':session_state['db']['dec'], 
                        'Name':session_state['db']['simbad_name'],
                        'RA':session_state['db']['ra_hms'],
                        'Vmag':session_state['db']['Vmag'],
                        'Imag':session_state['db']['Imag'],
                        'Gaiagmag':session_state['db']['phot_g_mean_mag'],
                        'Index':session_state['db']['Num']
                               })
data=ColumnDataSource(data=datadf)

# mapper = log_cmap(field_name='color', 
#                          palette=Viridis256[::-1],
#                          #palette=Turbo256[::-1],
#                          low=min(dist), high=max(dist),
#                         #low_color=Magma256[150], high_color=Magma256[200]
#                         )

tools = "hover, zoom_in, zoom_out, save, undo, redo, pan"
tooltips = [
        ('Num', '@Index'),
        ('Name', '@Name'),
        ('RA', '@RA'),
        ('Vmag', '@Vmag'),
        ('Imag', '@Imag'),
        ('Gaia g mag', '@Gaiagmag'),
    ]
p = figure(x_axis_label='RA', y_axis_label='DEC',
        background_fill_color='#222831', border_fill_color='#31363F',outline_line_color='#31363F',
        tools=tools, 
        tooltips=tooltips, toolbar_location="above", width=800, height=400)
p.yaxis.major_label_text_color = "#EEEEEE"
p.yaxis.axis_label_text_color = "#EEEEEE"
p.xaxis.major_label_text_color = "#EEEEEE"
p.xaxis.axis_label_text_color = "#EEEEEE"
p.grid.grid_line_color = '#EEEEEE'


p.circle('plotx','ploty', source=data, fill_alpha=0.6, color="#EEEEEE", size=20)
# color_bar = ColorBar(color_mapper=mapper['transform'], width=8, 
#                          location=(0,0), title="Distance",
#                         title_text_font_size = '12pt',
#                          major_label_text_font_size = '10pt',
#                          background_fill_color='#222831',major_label_text_color = "#EEEEEE",
#                          title_text_color = "#EEEEEE")
# p.add_layout(color_bar, 'right')
#st.text('Marker size corresponds to MS g magnitude; marker color corresponds to distance')
st.bokeh_chart(p, use_container_width=True)



############################ Form to update database:

'''## Update database

I'm currently working on a way for users to update the database in real time, but for now if you observe a star and 
find that it's a multiple, please email me at lapearce@umich.edu (or on slack) with the star(s) name and number of components
 '''




    
    

# ''' ## Update Database'''

with st.form("update"):
    left_co, cent_co,last_co = st.columns(3)
    with left_co:
        name = st.text_input(r"Simbad Name", key='name')
    with cent_co:
        N_sys = st.text_input(r"N components", key='n_sys')
    with last_co:
        vetted = st.text_input(r"Vetted", key='vetted')
    drop = st.checkbox("Drop from db?")
    submit = st.form_submit_button('Update')

# if submit:
#     sqlstring = "UPDATE aodb SET n_sys = "+N_sys+" WHERE simbad_name LIKE '%"+name+"%'"
#     querySQL(sqlstring)

from sqlalchemy import update
if submit:
    sqlstring = "UPDATE aodb SET n_sys = "+N_sys+" WHERE simbad_name LIKE '%"+name+"%'"
    with engine.connect() as conn:
        result = conn.execute(text(sqlstring))

    session_state['db'] = pd.DataFrame(result)

session_state['db']



