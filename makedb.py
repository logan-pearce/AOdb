import sqlite3
import pandas as pd
AOdb = pd.read_csv('Bright-AO-Stars.csv')
conn = sqlite3.connect('aodb.db')
AOdb.to_sql('aodb.db',conn,index=False, if_exists='replace')

# import sqlite3
# conn = sqlite3.connect('slsdb.db')
# slsdb.to_sql('slsdb.db',conn,index=False, if_exists='replace')