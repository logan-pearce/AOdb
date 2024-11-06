import streamlit as st

from pathlib import Path
import streamlit as st

# def read_markdown_file(markdown_file):
#     return Path(markdown_file).read_text()

# intro_markdown = read_markdown_file("Bright-Stars-Select.md")
# st.markdown(intro_markdown, unsafe_allow_html=True)
# How Bright Stars Were Selected

st.markdown("""I first performed a Gaia ADQL query to select bright stars within a dec range with RUWE < 1.2. This RUWE cut does not guarantee single stars since it only probes a subset of binary separation/luminosity ratio space (see El-Badry et al. 2024 Fig 1), but it will decrease the number likely to be on the edge of resolved and look blob-y on sky.

# Query
Query all objects in Gaia with continuous spectra, low RUWE, and within 15 deg of LCO declination. Include only one with continuous spectra to enable converying the Gaia spectra into WFS magnitudes.

I ran this query on the Gaia archive ADQL interface (because python TAP queries have a 2000 source limit). 


```python
query_input = '''SELECT source_id, ra, dec, pmra, pmdec, phot_g_mean_mag, bp_rp, ruwe, has_xp_continuous
                FROM gaiadr3.gaia_source WHERE
                has_xp_continuous = 'True' AND
                dec < -5 AND dec > -55 AND
                phot_g_mean_mag < 5 AND
                RUWE < 1.2
                '''

#from astroquery.gaia import Gaia
#job = Gaia.launch_job(query_input)
#r = job.get_results()

#### Load results from Gaia archive ADQL interface. 
r = pd.read_csv('1730739299568O-result.csv', dtype ={'source_id':str})
r
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>source_id</th>
      <th>ra</th>
      <th>dec</th>
      <th>pmra</th>
      <th>pmdec</th>
      <th>parallax</th>
      <th>phot_g_mean_mag</th>
      <th>bp_rp</th>
      <th>ruwe</th>
      <th>has_xp_continuous</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>5899408305748518528</td>
      <td>219.471578</td>
      <td>-49.425958</td>
      <td>-28.920706</td>
      <td>-29.213251</td>
      <td>9.348442</td>
      <td>4.001947</td>
      <td>-0.190720</td>
      <td>1.141235</td>
      <td>True</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2940472157174944128</td>
      <td>99.472556</td>
      <td>-18.237520</td>
      <td>-7.165916</td>
      <td>-8.777130</td>
      <td>6.787632</td>
      <td>4.103014</td>
      <td>1.312503</td>
      <td>1.031607</td>
      <td>True</td>
    </tr>
    <tr>
      <th>2</th>
      <td>5518235833732209536</td>
      <td>118.325618</td>
      <td>-48.102908</td>
      <td>-5.679107</td>
      <td>6.162858</td>
      <td>1.528783</td>
      <td>4.174099</td>
      <td>-0.184238</td>
      <td>1.048019</td>
      <td>True</td>
    </tr>
    <tr>
      <th>3</th>
      <td>5440991156144997504</td>
      <td>158.803364</td>
      <td>-39.562578</td>
      <td>-31.504189</td>
      <td>2.619350</td>
      <td>3.615400</td>
      <td>4.187924</td>
      <td>2.791997</td>
      <td>1.030684</td>
      <td>True</td>
    </tr>
    <tr>
      <th>4</th>
      <td>4825052919883631616</td>
      <td>76.102363</td>
      <td>-35.483166</td>
      <td>125.971691</td>
      <td>-42.909247</td>
      <td>17.591152</td>
      <td>4.194447</td>
      <td>1.369841</td>
      <td>1.109602</td>
      <td>True</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>188</th>
      <td>6816586167927675776</td>
      <td>325.503348</td>
      <td>-23.263280</td>
      <td>91.684025</td>
      <td>-95.142902</td>
      <td>13.485332</td>
      <td>4.984811</td>
      <td>1.128248</td>
      <td>0.939557</td>
      <td>True</td>
    </tr>
    <tr>
      <th>189</th>
      <td>4112977922408078848</td>
      <td>254.990472</td>
      <td>-25.092270</td>
      <td>5.211071</td>
      <td>-15.899578</td>
      <td>6.048476</td>
      <td>4.991735</td>
      <td>2.204937</td>
      <td>1.086619</td>
      <td>True</td>
    </tr>
    <tr>
      <th>190</th>
      <td>3011410242215060224</td>
      <td>89.768015</td>
      <td>-9.558458</td>
      <td>9.575414</td>
      <td>-46.119938</td>
      <td>11.299821</td>
      <td>4.998930</td>
      <td>0.238009</td>
      <td>1.071681</td>
      <td>True</td>
    </tr>
    <tr>
      <th>191</th>
      <td>6098305973471310080</td>
      <td>219.079095</td>
      <td>-46.245374</td>
      <td>-31.247339</td>
      <td>17.232811</td>
      <td>4.698888</td>
      <td>4.999646</td>
      <td>1.669282</td>
      <td>1.186192</td>
      <td>True</td>
    </tr>
    <tr>
      <th>192</th>
      <td>5416714660964273536</td>
      <td>153.881629</td>
      <td>-43.112611</td>
      <td>21.955339</td>
      <td>-54.174673</td>
      <td>6.037673</td>
      <td>4.999767</td>
      <td>1.735911</td>
      <td>1.050030</td>
      <td>True</td>
    </tr>
  </tbody>
</table>
<p>193 rows × 10 columns</p>
</div>



# Generate Synthetic Photometry using GaiaXpy

Next I used photometry to convert Gaia spectra into SDSS photometry to get WFS magnitude estimates.

Run the query on the Gaia archive, download the results and the continuous raw spectra.

Gaia expresses the BP/RP spectra as a function of coefficients of basis functions.  So to use the spectra you have to use the GaiaXPy tool to convert the continuous spectra into sampled spectra at a range of wavelengths.  

The convert function returns sampled spectra in pseudo-units (e s$^{-1}$ and pseudo-wavelengths)

The calibrate function returns sampled spectra in absolute units (W nm$^{-1}$ m$^{-2}$ and nm)

The generate function retrieves synthetic photometry in the filter set of your choice


https://gaia-dpci.github.io/GaiaXPy-website/tutorials.html

# Run the above query in the Gaia Archive ADQL interface
# Download the query results


![pic1](images/pic1.png)

# And the continuous spectra:

![pic2](images/pic2.png)

# And the sampled spectra:

![pic3](images/pic3.png)

### These files (and this notebook) are included in the Github repo for this app to enable reproducing this catalog.

# Use the continuous raw spectra to generate synthetic SDSS i band photometry:


```python
from gaiaxpy import generate, PhotometricSystem

# Path to file with XP CONTINUOUS RAW data (csv, ecsv, fits, or xml)
f = 'XP_CONTINUOUS_RAW.csv'
# Select a photometric system
phot_system = [PhotometricSystem.SDSS,PhotometricSystem.SDSS_Std]

# Generate snythetic phometry in SDSS system:
synthetic_photometry = generate(f, photometric_system=phot_system)

```

    Done! Output saved to path: ./output_synthetic_photometry.csv                   193 [00:00<?, ?spec/s][0m


```python
synthetic_photometry.to_csv('Bright-AO-Stars-Synthetic-Photometry.csv', index=False)
synthetic_photometry
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>source_id</th>
      <th>Sdss_mag_u</th>
      <th>Sdss_mag_g</th>
      <th>Sdss_mag_r</th>
      <th>Sdss_mag_i</th>
      <th>Sdss_mag_z</th>
      <th>Sdss_flux_u</th>
      <th>Sdss_flux_g</th>
      <th>Sdss_flux_r</th>
      <th>Sdss_flux_i</th>
      <th>...</th>
      <th>SdssStd_flux_u</th>
      <th>SdssStd_flux_g</th>
      <th>SdssStd_flux_r</th>
      <th>SdssStd_flux_i</th>
      <th>SdssStd_flux_z</th>
      <th>SdssStd_flux_error_u</th>
      <th>SdssStd_flux_error_g</th>
      <th>SdssStd_flux_error_r</th>
      <th>SdssStd_flux_error_i</th>
      <th>SdssStd_flux_error_z</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>5989102478619519616</td>
      <td>5.755067</td>
      <td>4.753345</td>
      <td>4.562233</td>
      <td>4.522348</td>
      <td>4.570907</td>
      <td>1.811228e-25</td>
      <td>4.556822e-25</td>
      <td>5.433838e-25</td>
      <td>5.637166e-25</td>
      <td>...</td>
      <td>2.381881e-25</td>
      <td>4.577204e-25</td>
      <td>5.428518e-25</td>
      <td>5.633938e-25</td>
      <td>5.435555e-25</td>
      <td>3.016513e-27</td>
      <td>4.138486e-28</td>
      <td>5.804783e-28</td>
      <td>2.520995e-28</td>
      <td>4.309222e-28</td>
    </tr>
    <tr>
      <th>1</th>
      <td>3507879565090229888</td>
      <td>6.377998</td>
      <td>5.030770</td>
      <td>4.529934</td>
      <td>4.376079</td>
      <td>4.327418</td>
      <td>1.020471e-25</td>
      <td>3.529329e-25</td>
      <td>5.597915e-25</td>
      <td>6.450128e-25</td>
      <td>...</td>
      <td>1.272854e-25</td>
      <td>3.585688e-25</td>
      <td>5.588167e-25</td>
      <td>6.459192e-25</td>
      <td>6.802406e-25</td>
      <td>1.141681e-27</td>
      <td>2.200916e-28</td>
      <td>6.117793e-28</td>
      <td>3.305058e-28</td>
      <td>6.124760e-28</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2944471321481716224</td>
      <td>7.978153</td>
      <td>5.797811</td>
      <td>4.747534</td>
      <td>4.387583</td>
      <td>4.171607</td>
      <td>2.337431e-26</td>
      <td>1.741308e-25</td>
      <td>4.581277e-25</td>
      <td>6.382148e-25</td>
      <td>...</td>
      <td>2.138519e-26</td>
      <td>1.805745e-25</td>
      <td>4.564387e-25</td>
      <td>6.413545e-25</td>
      <td>7.827252e-25</td>
      <td>3.558263e-28</td>
      <td>1.053251e-28</td>
      <td>3.628949e-28</td>
      <td>2.468492e-28</td>
      <td>4.788095e-28</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3185689878960742144</td>
      <td>8.366510</td>
      <td>6.021334</td>
      <td>4.741537</td>
      <td>4.257475</td>
      <td>3.962593</td>
      <td>1.634542e-26</td>
      <td>1.417315e-25</td>
      <td>4.606651e-25</td>
      <td>7.194653e-25</td>
      <td>...</td>
      <td>1.284329e-26</td>
      <td>1.483461e-25</td>
      <td>4.584611e-25</td>
      <td>7.243950e-25</td>
      <td>9.518014e-25</td>
      <td>2.924330e-28</td>
      <td>9.747137e-29</td>
      <td>3.790675e-28</td>
      <td>2.399517e-28</td>
      <td>5.528253e-28</td>
    </tr>
    <tr>
      <th>4</th>
      <td>6477431326619860224</td>
      <td>7.605403</td>
      <td>5.604990</td>
      <td>4.718381</td>
      <td>4.418542</td>
      <td>4.246637</td>
      <td>3.294873e-26</td>
      <td>2.079716e-25</td>
      <td>4.705953e-25</td>
      <td>6.202732e-25</td>
      <td>...</td>
      <td>3.326067e-26</td>
      <td>2.140902e-25</td>
      <td>4.691647e-25</td>
      <td>6.226528e-25</td>
      <td>7.312160e-25</td>
      <td>4.168134e-28</td>
      <td>1.103400e-28</td>
      <td>3.630606e-28</td>
      <td>2.160705e-28</td>
      <td>4.725068e-28</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>188</th>
      <td>2924859847973986048</td>
      <td>4.036712</td>
      <td>4.044914</td>
      <td>4.530514</td>
      <td>4.842717</td>
      <td>5.143209</td>
      <td>8.816883e-25</td>
      <td>8.750527e-25</td>
      <td>5.594928e-25</td>
      <td>4.196753e-25</td>
      <td>...</td>
      <td>1.452590e-24</td>
      <td>8.688678e-25</td>
      <td>5.599257e-25</td>
      <td>4.171597e-25</td>
      <td>3.193073e-25</td>
      <td>1.539558e-26</td>
      <td>9.418193e-28</td>
      <td>9.371379e-28</td>
      <td>1.573253e-28</td>
      <td>2.315395e-28</td>
    </tr>
    <tr>
      <th>189</th>
      <td>5699115357749939328</td>
      <td>7.512182</td>
      <td>5.628330</td>
      <td>4.761106</td>
      <td>4.485889</td>
      <td>4.329223</td>
      <td>3.590272e-26</td>
      <td>2.035485e-25</td>
      <td>4.524363e-25</td>
      <td>5.829675e-25</td>
      <td>...</td>
      <td>3.697489e-26</td>
      <td>2.093624e-25</td>
      <td>4.511274e-25</td>
      <td>5.848738e-25</td>
      <td>6.760658e-25</td>
      <td>3.437738e-28</td>
      <td>9.226239e-29</td>
      <td>2.990512e-28</td>
      <td>1.793591e-28</td>
      <td>3.962886e-28</td>
    </tr>
    <tr>
      <th>190</th>
      <td>4109030160308320128</td>
      <td>7.015432</td>
      <td>5.486960</td>
      <td>4.820432</td>
      <td>4.595452</td>
      <td>4.502999</td>
      <td>5.673186e-26</td>
      <td>2.318547e-25</td>
      <td>4.283779e-25</td>
      <td>5.270103e-25</td>
      <td>...</td>
      <td>6.515238e-26</td>
      <td>2.372309e-25</td>
      <td>4.272858e-25</td>
      <td>5.282430e-25</td>
      <td>5.784202e-25</td>
      <td>8.856953e-28</td>
      <td>1.979881e-28</td>
      <td>5.374413e-28</td>
      <td>2.672068e-28</td>
      <td>5.424166e-28</td>
    </tr>
    <tr>
      <th>191</th>
      <td>6521403820271516160</td>
      <td>7.639176</td>
      <td>5.686427</td>
      <td>4.786995</td>
      <td>4.470931</td>
      <td>4.289656</td>
      <td>3.193961e-26</td>
      <td>1.929432e-25</td>
      <td>4.417760e-25</td>
      <td>5.910547e-25</td>
      <td>...</td>
      <td>3.183425e-26</td>
      <td>1.988762e-25</td>
      <td>4.403142e-25</td>
      <td>5.934159e-25</td>
      <td>7.029398e-25</td>
      <td>3.367026e-28</td>
      <td>9.534291e-29</td>
      <td>3.198584e-28</td>
      <td>2.175558e-28</td>
      <td>4.237760e-28</td>
    </tr>
    <tr>
      <th>192</th>
      <td>5158225936898617344</td>
      <td>4.631533</td>
      <td>4.019388</td>
      <td>4.394299</td>
      <td>4.631626</td>
      <td>4.843334</td>
      <td>5.097848e-25</td>
      <td>8.958697e-25</td>
      <td>6.342792e-25</td>
      <td>5.097412e-25</td>
      <td>...</td>
      <td>7.730373e-25</td>
      <td>8.904086e-25</td>
      <td>6.345090e-25</td>
      <td>5.072511e-25</td>
      <td>4.225788e-25</td>
      <td>9.424578e-27</td>
      <td>9.109226e-28</td>
      <td>9.569257e-28</td>
      <td>2.920971e-28</td>
      <td>4.079733e-28</td>
    </tr>
  </tbody>
</table>
<p>193 rows × 31 columns</p>
</div>



# Using individual sampled spectra to generate magnitude in MagAO-X WFS

## Example calculation:

### Load MagAO-X WFS filter curves and Gaia G band curve:


```python
f6535 = pd.read_table('filter_curves/magaox_wfs-open_bs-65-35_atm.dat', comment='#', 
                  names=['wavelength [m]','transmission'], sep='\s+')
f6535['normalized transmission'] = f6535['transmission']/np.max(f6535['transmission'])
f6535['wavelength [nm]'] = f6535['wavelength [m]']*u.m.to(u.nm)

fhair = pd.read_table('filter_curves/magaox_wfs-open_bs-halpha-ir.dat', comment='#', 
                  names=['wavelength [m]','transmission'], sep='\s+')
fhair['normalized transmission'] = fhair['transmission']/np.max(fhair['transmission'])
fhair['wavelength [nm]'] = fhair['wavelength [m]']*u.m.to(u.nm)

g = pd.read_table('filter_curves/GaiaEDR3_passbands_zeropoints_version2/passband.dat', comment='#', 
                  names=['wavelength [nm]','G','eG','BP','eBP','RP','eRP'], sep='\s+')
g.loc[np.where(g['BP']==99.99)[0], 'BP'] = 0
g.loc[np.where(g['RP']==99.99)[0], 'RP'] = 0
g.loc[np.where(g['G']==99.99)[0], 'G'] = 0
g['normalized transmission'] = g['G']/np.max(g['G'])

%matplotlib inline
plt.plot(g['wavelength [nm]'], g['normalized transmission'], label='Gaia g',color='darkviolet')
plt.plot(f6535['wavelength [nm]'],f6535['normalized transmission'], label='WFS 65/35',color='blue')
plt.plot(fhair['wavelength [nm]'],fhair['normalized transmission'], 
         label=r'WFS H$\alpha$/IR',color='orange')
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Transmission')
plt.legend(fontsize=15,loc='upper left')
plt.tight_layout()
plt.grid(ls=':')

```


    
![png](images/output_9_0.png)
    


### Load the Gaia sampled spectrum downloaded from the archive for a single source:


```python
# load individual spectrum:
single_source_spectrum = pd.read_csv('lpearce1730735985448-sampled-spectra/XP_SAMPLED-Gaia DR3 3622819624439843072.csv',
                                    dtype={'source_id':str})
single_source_spectrum
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>source_id</th>
      <th>solution_id</th>
      <th>ra</th>
      <th>dec</th>
      <th>wavelength</th>
      <th>flux</th>
      <th>flux_error</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>3622819624439843072</td>
      <td>4545469030156206080</td>
      <td>196.974296</td>
      <td>-10.740454</td>
      <td>336.0</td>
      <td>8.989019e-14</td>
      <td>1.449126e-14</td>
    </tr>
    <tr>
      <th>1</th>
      <td>3622819624439843072</td>
      <td>4545469030156206080</td>
      <td>196.974296</td>
      <td>-10.740454</td>
      <td>338.0</td>
      <td>8.974217e-14</td>
      <td>1.095063e-14</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3622819624439843072</td>
      <td>4545469030156206080</td>
      <td>196.974296</td>
      <td>-10.740454</td>
      <td>340.0</td>
      <td>8.739168e-14</td>
      <td>8.897407e-15</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3622819624439843072</td>
      <td>4545469030156206080</td>
      <td>196.974296</td>
      <td>-10.740454</td>
      <td>342.0</td>
      <td>7.167989e-14</td>
      <td>7.729005e-15</td>
    </tr>
    <tr>
      <th>4</th>
      <td>3622819624439843072</td>
      <td>4545469030156206080</td>
      <td>196.974296</td>
      <td>-10.740454</td>
      <td>344.0</td>
      <td>5.663622e-14</td>
      <td>7.209097e-15</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>338</th>
      <td>3622819624439843072</td>
      <td>4545469030156206080</td>
      <td>196.974296</td>
      <td>-10.740454</td>
      <td>1012.0</td>
      <td>2.126632e-13</td>
      <td>9.147954e-15</td>
    </tr>
    <tr>
      <th>339</th>
      <td>3622819624439843072</td>
      <td>4545469030156206080</td>
      <td>196.974296</td>
      <td>-10.740454</td>
      <td>1014.0</td>
      <td>2.085074e-13</td>
      <td>1.029905e-14</td>
    </tr>
    <tr>
      <th>340</th>
      <td>3622819624439843072</td>
      <td>4545469030156206080</td>
      <td>196.974296</td>
      <td>-10.740454</td>
      <td>1016.0</td>
      <td>2.122904e-13</td>
      <td>1.115871e-14</td>
    </tr>
    <tr>
      <th>341</th>
      <td>3622819624439843072</td>
      <td>4545469030156206080</td>
      <td>196.974296</td>
      <td>-10.740454</td>
      <td>1018.0</td>
      <td>2.122881e-13</td>
      <td>1.096124e-14</td>
    </tr>
    <tr>
      <th>342</th>
      <td>3622819624439843072</td>
      <td>4545469030156206080</td>
      <td>196.974296</td>
      <td>-10.740454</td>
      <td>1020.0</td>
      <td>2.213856e-13</td>
      <td>1.051290e-14</td>
    </tr>
  </tbody>
</table>
<p>343 rows × 7 columns</p>
</div>



### Find which source this is in the catalog of downloaded sources:


```python
ind = np.where(r['source_id'] == single_source_spectrum['source_id'][0])[0][0]
ind
```




    137



### Compute $\lambda_0$ for MagAO-X WFS and Gaia G filter:

$$\lambda_0 = \frac{\int_{0}^{\infty} \lambda R(\lambda) d\lambda}{\int_{0}^{\infty} R(\lambda) d\lambda}
$$ 
Where $R(\lambda)$ is the filter response curve


```python
# Determine step sizes for each filter set:
dl_f6535 = np.mean([f6535['wavelength [nm]'][i+1] - f6535['wavelength [nm]'][i] 
           for i in range(1,len(f6535['wavelength [nm]'])-1)])
dl_fhair = np.mean([fhair['wavelength [nm]'][i+1] - fhair['wavelength [nm]'][i] 
           for i in range(1,len(fhair['wavelength [nm]'])-1)])
dl_g = np.mean([g['wavelength [nm]'][i+1] - g['wavelength [nm]'][i] 
           for i in range(1,len(g['wavelength [nm]'])-1)])

```


```python
GaiaG_lambda0 = np.sum(g['wavelength [nm]']*g['normalized transmission']* dl_g)/ \
            np.sum(g['normalized transmission']* dl_g)
f6535_lambda0 = np.sum(f6535['wavelength [nm]']*f6535['normalized transmission']* dl_f6535)/ \
            np.sum(f6535['normalized transmission']* dl_f6535)
fhair_lambda0 = np.sum(fhair['wavelength [nm]']*fhair['normalized transmission']* dl_fhair)/ \
            np.sum(fhair['normalized transmission']* dl_fhair)
```

## Compute flux at effective wavelength:

The flux at $\lambda_0$ is the integral of the flux times the filter response curve divided by the effective width of the filter

$$ F_\lambda(\lambda_0) \Delta\lambda = \int_{0}^{\infty} F_\lambda(\lambda) R(\lambda) d\lambda $$


```python
from scipy.interpolate import interp1d

def FluxLambda0(spectrum_wavelength,spectrum_flux, lambda_0):
    from scipy.interpolate import interp1d
    # create the interpolation function:
    interpfunc = interp1d(spectrum_wavelength,spectrum_flux, fill_value="extrapolate")
    # Interpolate the filter's central wavelength in the spectrum's flux array:
    F_lambda_0 = interpfunc(lambda_0)
    return F_lambda_0

# Get the flux at Lambda_0
F_lambda0_g = FluxLambda0(single_source_spectrum['wavelength'], single_source_spectrum['flux'], GaiaG_lambda0)
F_lambda0_f6535 = FluxLambda0(single_source_spectrum['wavelength'], single_source_spectrum['flux'], f6535_lambda0)
F_lambda0_fhair = FluxLambda0(single_source_spectrum['wavelength'], single_source_spectrum['flux'], fhair_lambda0)

# Compute for Vega:
vega = pd.read_csv('vega.csv')
vega_F_lambda0_g = FluxLambda0(vega['WAVELENGTH']*u.AA.to(u.nm), vega['FLUX'], GaiaG_lambda0)
vega_F_lambda0_f6535 = FluxLambda0(vega['WAVELENGTH']*u.AA.to(u.nm), vega['FLUX'], f6535_lambda0)
vega_F_lambda0_fhair = FluxLambda0(vega['WAVELENGTH']*u.AA.to(u.nm), vega['FLUX'], fhair_lambda0)
```


```python
# Colors:
GaiaG_to_f6535 = -2.5*np.log10(F_lambda0_g/vega_F_lambda0_g) - (-2.5*np.log10(F_lambda0_f6535/vega_F_lambda0_f6535))
GaiaG_to_fhair = -2.5*np.log10(F_lambda0_g/vega_F_lambda0_g) - (-2.5*np.log10(F_lambda0_fhair/vega_F_lambda0_fhair))

print('G_to_f6535',GaiaG_to_f6535,'G_to_fhair',GaiaG_to_fhair)
print('Mag in G:',r['phot_g_mean_mag'][ind])
print('Mag in 6535 WFS:',r['phot_g_mean_mag'][ind]-GaiaG_to_f6535)
print('Mag in HaIR WFS:',r['phot_g_mean_mag'][ind]-GaiaG_to_fhair)
```

    G_to_f6535 0.5775996185317158 G_to_fhair 0.6840397501937421
    Mag in G: 4.8242455
    Mag in 6535 WFS: 4.246645881468284
    Mag in HaIR WFS: 4.140205749806258


# Do it for all the sources:


```python
r['WFS6535 mag'] = np.nan
r['WFSHaIR mag'] = np.nan
```


```python
def FluxLambda0(spectrum_wavelength,spectrum_flux, lambda_0):
    from scipy.interpolate import interp1d
    # create the interpolation function:
    interpfunc = interp1d(spectrum_wavelength,spectrum_flux, fill_value="extrapolate")
    # Interpolate the filter's central wavelength in the spectrum's flux array:
    F_lambda_0 = interpfunc(lambda_0)
    return F_lambda_0

os.system('ls lpearce1730735985448-sampled-spectra/* > list')
with open('list') as f:
    z = f.read().splitlines()

f6535 = pd.read_table('filter_curves/magaox_wfs-open_bs-65-35_atm.dat', comment='#', 
                  names=['wavelength [m]','transmission'], sep='\s+')
f6535['normalized transmission'] = f6535['transmission']/np.max(f6535['transmission'])
f6535['wavelength [nm]'] = f6535['wavelength [m]']*u.m.to(u.nm)

fhair = pd.read_table('filter_curves/magaox_wfs-open_bs-halpha-ir.dat', comment='#', 
                  names=['wavelength [m]','transmission'], sep='\s+')
fhair['normalized transmission'] = fhair['transmission']/np.max(fhair['transmission'])
fhair['wavelength [nm]'] = fhair['wavelength [m]']*u.m.to(u.nm)

g = pd.read_table('filter_curves/GaiaEDR3_passbands_zeropoints_version2/passband.dat', comment='#', 
                  names=['wavelength [nm]','G','eG','BP','eBP','RP','eRP'], sep='\s+')
g.loc[np.where(g['BP']==99.99)[0], 'BP'] = 0
g.loc[np.where(g['RP']==99.99)[0], 'RP'] = 0
g.loc[np.where(g['G']==99.99)[0], 'G'] = 0
g['normalized transmission'] = g['G']/np.max(g['G'])
```


```python
k = pd.read_csv(z[0])
# Determine step sizes for each filter set:
dl_f6535 = np.mean([f6535['wavelength [nm]'][i+1] - f6535['wavelength [nm]'][i] 
           for i in range(1,len(f6535['wavelength [nm]'])-1)])
dl_fhair = np.mean([fhair['wavelength [nm]'][i+1] - fhair['wavelength [nm]'][i] 
           for i in range(1,len(fhair['wavelength [nm]'])-1)])
dl_g = np.mean([g['wavelength [nm]'][i+1] - g['wavelength [nm]'][i] 
           for i in range(1,len(g['wavelength [nm]'])-1)])

# Determine effective wavelength:
GaiaG_lambda0 = np.sum(g['wavelength [nm]']*g['normalized transmission']* dl_g)/ \
            np.sum(g['normalized transmission']* dl_g)
f6535_lambda0 = np.sum(f6535['wavelength [nm]']*f6535['normalized transmission']* dl_f6535)/ \
            np.sum(f6535['normalized transmission']* dl_f6535)
fhair_lambda0 = np.sum(fhair['wavelength [nm]']*fhair['normalized transmission']* dl_fhair)/ \
            np.sum(fhair['normalized transmission']* dl_fhair)

# Compute for Vega:
vega = pd.read_csv('vega.csv')
vega_F_lambda0_g = FluxLambda0(vega['WAVELENGTH']*u.AA.to(u.nm), vega['FLUX'], GaiaG_lambda0)
vega_F_lambda0_f6535 = FluxLambda0(vega['WAVELENGTH']*u.AA.to(u.nm), vega['FLUX'], f6535_lambda0)
vega_F_lambda0_fhair = FluxLambda0(vega['WAVELENGTH']*u.AA.to(u.nm), vega['FLUX'], fhair_lambda0)
```


```python
from scipy.interpolate import interp1d
from myastrotools.tools import update_progress

for j,spectrum in enumerate(z):
    k = pd.read_csv(spectrum, dtype={'source_id':str})
    i = np.where(r['source_id'] == k['source_id'][0])[0][0]
    
    # make interpolation function:
    spl = interp1d(k['wavelength'], k['flux'])
    
    # interpolate:
    F_lambda0_g = spl(GaiaG_lambda0)
    F_lambda0_f6535 = spl(f6535_lambda0)
    F_lambda0_fhair = spl(fhair_lambda0)
    
    # compute colors:
    GaiaG_to_f6535 = -2.5*np.log10(F_lambda0_g/vega_F_lambda0_g) - (-2.5*np.log10(F_lambda0_f6535/vega_F_lambda0_f6535))
    GaiaG_to_fhair = -2.5*np.log10(F_lambda0_g/vega_F_lambda0_g) - (-2.5*np.log10(F_lambda0_fhair/vega_F_lambda0_fhair))
    
    # compute mags in filters:
    r.loc[i, 'WFS6535 mag'] = r['phot_g_mean_mag'][i]-GaiaG_to_f6535
    r.loc[i, 'WFSHaIR mag'] = r['phot_g_mean_mag'][i]-GaiaG_to_fhair
    
    update_progress(j,len(z))
```

    99.0% (192 of 193): |####################|  

### Put the GaiaXPy mags in the table:


```python
synthetic_photometry = pd.read_csv('Bright-AO-Stars-Synthetic-Photometry.csv', dtype={'source_id':str})
cols = synthetic_photometry.columns
cols
```




    Index(['source_id', 'Sdss_mag_u', 'Sdss_mag_g', 'Sdss_mag_r', 'Sdss_mag_i',
           'Sdss_mag_z', 'Sdss_flux_u', 'Sdss_flux_g', 'Sdss_flux_r',
           'Sdss_flux_i', 'Sdss_flux_z', 'Sdss_flux_error_u', 'Sdss_flux_error_g',
           'Sdss_flux_error_r', 'Sdss_flux_error_i', 'Sdss_flux_error_z',
           'SdssStd_mag_u', 'SdssStd_mag_g', 'SdssStd_mag_r', 'SdssStd_mag_i',
           'SdssStd_mag_z', 'SdssStd_flux_u', 'SdssStd_flux_g', 'SdssStd_flux_r',
           'SdssStd_flux_i', 'SdssStd_flux_z', 'SdssStd_flux_error_u',
           'SdssStd_flux_error_g', 'SdssStd_flux_error_r', 'SdssStd_flux_error_i',
           'SdssStd_flux_error_z'],
          dtype='object')




```python
for col in cols[1:]:
    r[col] = np.nan

for j in range(len(r)):
    i = np.where(r['source_id'] == synthetic_photometry['source_id'][j])[0][0]
    for col in cols:
        r.loc[i, col] = synthetic_photometry[col][j]
        r.loc[i, col] = synthetic_photometry[col][j]
```

## Identify multiples from CMD

As another check for eliminating multiples, I used the Gai color-magnitude diagram to help select overluminous sources.


```python
def GetPointsWithinARegion(xdata, ydata, points):
    ''' For a region defined by points, return the indicies of items from [xdata,ydata]
    that lie within that region
    
    Args:
        xdata, ydata (arr): x and y data 
        points (arr): array of points describing region in tuples of (x,y)
        
    Returns:
        indicies of points in dataframe that lie within the region.
    '''
    y = points[:,1]
    x = points[:,0]

    # find points that lie within region:
    stacked1 = np.stack((xdata,ydata),axis=1)
    from matplotlib import path
    p = path.Path(points)
    indicieswithinregion = p.contains_points(stacked1)
    return indicieswithinregion
```

For comparison, I used a CMG of stars in the Praesepe cluster from Deacon & Kraus 2021 Fig 4 (http://adsabs.harvard.edu/abs/2020MNRAS.496.5176D). This is Gaia cmd of known members of the young stellar association (and thus almost represents a ZAMS) with high RUWE and excess noise sources (indicating multiplicity) identified with triangles.  

Below I plotted the Praesepe members in orange. Multiple systems lie above the single star main sequence since they are brighter than a single star. This is born out by the high RUQE and excess noise source lying mostly above the MS.  

I plotted stars in the catalog in blue circles. I selected catalog sources within the teal region and removed them from the catalog.  This is not a fool proof method - I drew the region by eye, and likey missed many multiples and possibly discarded some singles; it also does not address possible multiples in the RGB. But it is used as a rough pass for eliminating many multiples quickly.


```python
cmd = pickle.load(open('praesepe_cmd_data.pkl','rb'))
bprp, GMag, ol, an = cmd[0], cmd[1], cmd[2], cmd[3]
```


```python
%matplotlib inline
fig, ax = plt.subplots(figsize=(12,8))
ax.plot(bprp, GMag, marker=".", ls='None',markersize=1.5, color='orange')
ax.scatter(bprp[ol], GMag[ol], marker="^", ls='None',s=15, facecolors='None',edgecolors='blue')
ax.scatter(bprp[an], GMag[an], marker="v", ls='None',s=15, facecolors='None',edgecolors='purple')


ax.scatter(r['bp_rp'], 
          r['phot_g_mean_mag'] - 5*np.log10(1000/r['parallax']) + 5,
          alpha = 0.5)



ax.invert_yaxis()
ax.set_xlim(-1,4)
ax.set_xlabel('BP-RP [mag]')
ax.set_ylabel('Abs G Mag')
ax.grid(ls=':')
#ax.legend(loc='lower left', fontsize=20)

ax = plt.gca()
axins = ax.inset_axes([0.9, 0.3, 0.8, 0.8])

axins.plot(bprp, GMag, marker=".", ls='None',markersize=5, color='orange')
axins.scatter(bprp[ol], GMag[ol], marker="^", ls='None',s=40, facecolors='None',edgecolors='blue',zorder=10)
axins.scatter(bprp[an], GMag[an], marker="v", ls='None',s=40, facecolors='None',edgecolors='purple',zorder=10)
axins.scatter(r['bp_rp'], 
          r['phot_g_mean_mag'] - 5*np.log10(1000/r['parallax']) + 5,
          alpha = 0.5, zorder=10)


points = np.array([
    [0.85,4.5] , [0.15,0], [0, -3.5], [0.8,-4], [0.95, 4], [0.85,4.5]
])
multiples = GetPointsWithinARegion(r['bp_rp'], 
          r['phot_g_mean_mag'] - 5*np.log10(1000/r['parallax']) + 5, points)
axins.plot(points[:,0], points[:,1],color='teal')

axins.set_xlim(-0.5,1.25)
axins.set_ylim(-5,7)
axins.invert_yaxis()
axins.grid(ls=':')

plt.tight_layout()
#plt.savefig('praesepe_CMD.png',dpi=300,bbox_inches='tight')
```


    
![png](images/output_32_0.png)
    



```python
r_singles = r.loc[~multiples]
r_singles = r_singles.reset_index(drop=True)

```

Now I need to get more info about each source

### Get Simbad name, otype, and Vmag:


```python
# get simbad name:
from myastrotools.tools import update_progress
from astroquery.simbad import Simbad
import warnings
warnings.filterwarnings('ignore')
simbad = Simbad()
simbad.add_votable_fields('flux(V)','ids','sptype','otype')

r_singles['simbad_name'] = np.nan
r_singles['HD'] = np.nan
r_singles['Vmag'] = np.nan
r_singles['otype'] = np.nan
for i in range(len(r_singles)):
    try:
        o = 'Gaia DR3 '+str(r_singles['source_id'][i])
        obj = simbad.query_objects([o])
        otype = obj['OTYPE']
        ids = obj['IDS'][0].split('|')
        try:
            r_singles.loc[i,'HD'] = str(list(filter(lambda x: 'HD ' in x, ids))[0])
        except:
            pass
        r_singles.loc[i,'simbad_name'] = str(ids[0]).replace('* ','')
        r_singles.loc[i,'Vmag'] = obj['FLUX_V']
        r_singles.loc[i,'otype'] = otype

    except:
        pass

    update_progress(i,len(r_singles))
```

    99.0% (182 of 183): |####################|  

The 'otype' Simbad entry catalogs if a source is a known binary, so we can just drop those.


```python
mult = []
for i in range(len(r_singles)):
    if '**' in str(r_singles.loc[i,'otype']):
        mult.append(i)
r_singles = r_singles.drop(index=mult)
r_singles = r_singles.reset_index(drop=True)
```

### Cut out ones fainter than 5.5 mags in either WFS filter:


```python
r_final = r_singles.loc[np.where((r_singles['WFS6535 mag'] < 5.5) & (r_singles['WFSHaIR mag'] < 5.5))]
r_final = r_final.sort_values('ra')
r_final = r_final.reset_index(drop=True)
r_final.to_csv('Bright-AO-Stars.csv', index=False)
r_final
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>source_id</th>
      <th>ra</th>
      <th>dec</th>
      <th>pmra</th>
      <th>pmdec</th>
      <th>parallax</th>
      <th>phot_g_mean_mag</th>
      <th>bp_rp</th>
      <th>ruwe</th>
      <th>has_xp_continuous</th>
      <th>...</th>
      <th>SdssStd_flux_error_u</th>
      <th>SdssStd_flux_error_g</th>
      <th>SdssStd_flux_error_r</th>
      <th>SdssStd_flux_error_i</th>
      <th>SdssStd_flux_error_z</th>
      <th>simbad_name</th>
      <th>HD</th>
      <th>Vmag</th>
      <th>otype</th>
      <th>OTYPE</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>2479813756210251136</td>
      <td>26.496870</td>
      <td>-5.733429</td>
      <td>-10.221720</td>
      <td>-29.187017</td>
      <td>5.555357</td>
      <td>4.780030</td>
      <td>1.744846</td>
      <td>1.010017</td>
      <td>True</td>
      <td>...</td>
      <td>3.050719e-28</td>
      <td>1.184490e-28</td>
      <td>4.668909e-28</td>
      <td>4.034298e-28</td>
      <td>7.767560e-28</td>
      <td>HD  10824</td>
      <td>HD  10824</td>
      <td>5.340</td>
      <td>Star</td>
      <td>Star</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2461043306017773696</td>
      <td>27.395593</td>
      <td>-10.686824</td>
      <td>-147.867214</td>
      <td>-93.479828</td>
      <td>42.584807</td>
      <td>4.568229</td>
      <td>0.459817</td>
      <td>0.860792</td>
      <td>True</td>
      <td>...</td>
      <td>2.253731e-27</td>
      <td>3.909967e-28</td>
      <td>5.926557e-28</td>
      <td>2.484584e-28</td>
      <td>4.664472e-28</td>
      <td>chi Cet</td>
      <td>HD  11171</td>
      <td>4.680</td>
      <td>HighPM*</td>
      <td>HighPM*</td>
    </tr>
    <tr>
      <th>2</th>
      <td>5134765416778736512</td>
      <td>29.167801</td>
      <td>-22.526894</td>
      <td>60.374169</td>
      <td>-24.634717</td>
      <td>7.328994</td>
      <td>4.407405</td>
      <td>1.601481</td>
      <td>1.016402</td>
      <td>True</td>
      <td>...</td>
      <td>5.273653e-28</td>
      <td>1.932056e-28</td>
      <td>6.989209e-28</td>
      <td>5.122739e-28</td>
      <td>1.064459e-27</td>
      <td>56 Cet</td>
      <td>HD  11930</td>
      <td>4.850</td>
      <td>HighPM*</td>
      <td>HighPM*</td>
    </tr>
    <tr>
      <th>3</th>
      <td>5136714576016018432</td>
      <td>29.942590</td>
      <td>-20.824468</td>
      <td>18.826217</td>
      <td>15.859864</td>
      <td>4.763460</td>
      <td>4.662973</td>
      <td>2.047709</td>
      <td>0.835636</td>
      <td>True</td>
      <td>...</td>
      <td>3.173615e-28</td>
      <td>1.628789e-28</td>
      <td>6.022030e-28</td>
      <td>6.527044e-28</td>
      <td>1.199097e-27</td>
      <td>57 Cet</td>
      <td>HD  12255</td>
      <td>5.429</td>
      <td>Variable*</td>
      <td>Variable*</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5021010046848175616</td>
      <td>31.122699</td>
      <td>-29.296752</td>
      <td>1.849292</td>
      <td>11.772053</td>
      <td>9.403377</td>
      <td>4.641025</td>
      <td>-0.180284</td>
      <td>0.866612</td>
      <td>True</td>
      <td>...</td>
      <td>7.023842e-27</td>
      <td>6.832866e-28</td>
      <td>6.792276e-28</td>
      <td>1.581843e-28</td>
      <td>2.416016e-28</td>
      <td>HD  12767</td>
      <td>HD  12767</td>
      <td>4.690</td>
      <td>alf2CVnV*</td>
      <td>alf2CVnV*</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>174</th>
      <td>6628587855877381120</td>
      <td>335.398148</td>
      <td>-21.598607</td>
      <td>-10.986031</td>
      <td>-84.833784</td>
      <td>17.347473</td>
      <td>4.810523</td>
      <td>1.235893</td>
      <td>0.761455</td>
      <td>True</td>
      <td>...</td>
      <td>5.332644e-28</td>
      <td>1.669396e-28</td>
      <td>5.081345e-28</td>
      <td>3.144490e-28</td>
      <td>6.725791e-28</td>
      <td>47 Aqr</td>
      <td>HD 212010</td>
      <td>5.124</td>
      <td>RGB*</td>
      <td>RGB*</td>
    </tr>
    <tr>
      <th>175</th>
      <td>6601750220152445440</td>
      <td>337.876686</td>
      <td>-32.346159</td>
      <td>58.639817</td>
      <td>-18.857117</td>
      <td>22.083617</td>
      <td>4.263408</td>
      <td>0.027941</td>
      <td>1.097547</td>
      <td>True</td>
      <td>...</td>
      <td>4.844582e-27</td>
      <td>8.518736e-28</td>
      <td>8.559480e-28</td>
      <td>2.914520e-28</td>
      <td>4.536980e-28</td>
      <td>bet PsA</td>
      <td>HD 213398</td>
      <td>4.290</td>
      <td>HighPM*</td>
      <td>HighPM*</td>
    </tr>
    <tr>
      <th>176</th>
      <td>6544701440869855744</td>
      <td>340.874958</td>
      <td>-41.414762</td>
      <td>9.589843</td>
      <td>-93.616195</td>
      <td>13.415835</td>
      <td>4.550108</td>
      <td>1.187648</td>
      <td>1.020501</td>
      <td>True</td>
      <td>...</td>
      <td>4.917367e-28</td>
      <td>1.282133e-28</td>
      <td>3.597556e-28</td>
      <td>2.088870e-28</td>
      <td>4.053137e-28</td>
      <td>rho Gru</td>
      <td>HD 215104</td>
      <td>4.835</td>
      <td>HighPM*</td>
      <td>HighPM*</td>
    </tr>
    <tr>
      <th>177</th>
      <td>6546508900546867840</td>
      <td>342.759098</td>
      <td>-39.156869</td>
      <td>17.775482</td>
      <td>-7.583302</td>
      <td>3.627289</td>
      <td>4.930226</td>
      <td>1.580908</td>
      <td>1.008365</td>
      <td>True</td>
      <td>...</td>
      <td>3.695744e-28</td>
      <td>1.376046e-28</td>
      <td>4.210356e-28</td>
      <td>3.214627e-28</td>
      <td>6.583159e-28</td>
      <td>HD 216149</td>
      <td>HD 216149</td>
      <td>5.400</td>
      <td>Star</td>
      <td>Star</td>
    </tr>
    <tr>
      <th>178</th>
      <td>6521403820271516160</td>
      <td>359.732844</td>
      <td>-52.745531</td>
      <td>58.509812</td>
      <td>62.021589</td>
      <td>12.097462</td>
      <td>4.795882</td>
      <td>1.299198</td>
      <td>1.085166</td>
      <td>True</td>
      <td>...</td>
      <td>3.367026e-28</td>
      <td>9.534291e-29</td>
      <td>3.198584e-28</td>
      <td>2.175558e-28</td>
      <td>4.237760e-28</td>
      <td>HD 224554</td>
      <td>HD 224554</td>
      <td>5.133</td>
      <td>HighPM*</td>
      <td>HighPM*</td>
    </tr>
  </tbody>
</table>
<p>179 rows × 47 columns</p>
</div>



# Use Simbad for brighter stars not in Gaia

Gaia has an lower magnitude limit of ~3, so the brightest stars won't be in Gaia.  We can query Simbad directly for those.  We won't be able to estimate the WFS magnitudes though since we won't have a spectrum, but I band is close.  


```python
# Simbad TAP query:
from astroquery.simbad import Simbad
simbad = Simbad()
# Set max row limit to infinity so we get all the results:
simbad.ROW_LIMIT = -1
# The 'basic' table contains basic info on each star, the 'flux' table is the one with the magnitudes,
# so we can join them. Annoyingly, Simbad says "flux" when it means magnitudes.
# We can also automatically exclude binaries from the results (otype != **)
s = simbad.query_tap('''SELECT basic.ra, basic.dec, basic.oid, basic.pmra, basic.pmdec, main_id, flux, filter 
                FROM basic JOIN flux ON basic.oid = flux.oidref
                 WHERE (otype != '**') AND (flux < 3) AND (flux != 0.0)
                 AND (basic.dec < -5) AND (basic.dec > -55)
                 ''')
s = s.to_pandas()

```


```python
# How many have V mags, I mags, and i mags?
len(np.where(s['filter'] == 'V')[0]), len(np.where(s['filter'] == 'i')[0]), len(np.where(s['filter'] == 'I')[0]),
```




    (71, 3, 161)




```python
# Get V mags:
V = s.loc[np.where(s['filter'] == 'V')[0]]
V = V.reset_index(drop=True)
# V = V.drop([34,35])
# V = V.reset_index(drop=True)
V = V.sort_values('ra')
V = V.reset_index(drop=True)
```


```python
# Get I mags:
I = s.loc[np.where(s['filter'] == 'I')[0]]
I = I.sort_values('ra')
I = I.reset_index(drop=True)
```


```python
# Make a pandas dataframe of those:
brighter = pd.concat([V,I])
brighter = brighter.sort_values('ra')
brighter = brighter.reset_index(drop=True)
# Combine them into one.
brighter['Vmag'] = np.nan
brighter['Imag'] = np.nan
for i in range(len(brighter)):
    if brighter.loc[i,'filter'] == 'I':
        brighter.loc[i,'Imag'] = brighter.loc[i,'flux']
    elif brighter.loc[i,'filter'] == 'V':
        brighter.loc[i,'Vmag'] = brighter.loc[i,'flux']
len(brighter)
```




    232




```python
# Look for duplicate entries and drop one:
dupes = []
for i in range(len(brighter)):
    if len(np.where(brighter['oid'] == brighter['oid'][i])[0]) > 1:
        dupes.append(np.where(brighter['oid'] == brighter['oid'][i])[0])
        #print(np.where(brighter['oid'] == brighter['oid'][i])[0])
```


```python
for d in dupes[::2]:
    if np.isnan(brighter.loc[d[0]]['Vmag']):
        brighter.loc[d[0], 'Vmag'] = brighter.loc[d[1], 'Vmag']
    if np.isnan(brighter.loc[d[0]]['Imag']):
        brighter.loc[d[0], 'Imag'] = brighter.loc[d[1], 'Imag']
    brighter = brighter.drop(index=d[1])
brighter = brighter.reset_index(drop=True)
```


```python
# Double check we got them all:
dupes = []
for i in range(len(brighter)):
    if len(np.where(brighter['oid'] == brighter['oid'][i])[0]) > 1:
        dupes.append(np.where(brighter['oid'] == brighter['oid'][i])[0])
        #print(np.where(brighter['oid'] == brighter['oid'][i])[0])
dupes
```




    []




```python
# Simplfy some of the simbad names by removing the leading '*', which make future querying difficult.
for i in range(len(brighter)):
    brighter.loc[i,'main_id'] = str(brighter.loc[i,'main_id']).replace('* ','')
```


```python
# Eliminate the ones that don't have proper motions, since that's needed for the TCS catalog:
no_pm = []
for i in range(len(brighter)):
    if np.isnan(brighter.loc[i,'pmra']):
        no_pm.append(i)
    elif np.isnan(brighter.loc[i,'pmdec']):
        no_pm.append(i)

brighter = brighter.drop(index=no_pm)
brighter = brighter.reset_index(drop=True)
```

## Make useful database:


```python
r_out = r_final[
    ['simbad_name', 'ra', 'dec', 'Vmag', 'Sdss_mag_i', 'WFS6535 mag',
       'WFSHaIR mag', 'phot_g_mean_mag', 'pmra', 'pmdec', 'parallax',
       'bp_rp','source_id']
]
r_out['vetted'] = 'N'
r_out['n_sys'] = 1
r_out['Imag'] = r_out['Sdss_mag_i']
r_out = r_out[
    ['simbad_name', 'ra', 'dec', 'Vmag', 'Imag', 'WFS6535 mag',
       'WFSHaIR mag', 'phot_g_mean_mag', 'pmra', 'pmdec', 'parallax',
       'bp_rp','source_id']
]
r_out
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>simbad_name</th>
      <th>ra</th>
      <th>dec</th>
      <th>Vmag</th>
      <th>Imag</th>
      <th>WFS6535 mag</th>
      <th>WFSHaIR mag</th>
      <th>phot_g_mean_mag</th>
      <th>pmra</th>
      <th>pmdec</th>
      <th>parallax</th>
      <th>bp_rp</th>
      <th>source_id</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>HD  10824</td>
      <td>26.496870</td>
      <td>-5.733429</td>
      <td>5.340</td>
      <td>4.305772</td>
      <td>3.953154</td>
      <td>3.812200</td>
      <td>4.780030</td>
      <td>-10.221720</td>
      <td>-29.187017</td>
      <td>5.555357</td>
      <td>1.744846</td>
      <td>2479813756210251136</td>
    </tr>
    <tr>
      <th>1</th>
      <td>chi Cet</td>
      <td>27.395593</td>
      <td>-10.686824</td>
      <td>4.680</td>
      <td>4.632586</td>
      <td>4.325820</td>
      <td>4.295802</td>
      <td>4.568229</td>
      <td>-147.867214</td>
      <td>-93.479828</td>
      <td>42.584807</td>
      <td>0.459817</td>
      <td>2461043306017773696</td>
    </tr>
    <tr>
      <th>2</th>
      <td>56 Cet</td>
      <td>29.167801</td>
      <td>-22.526894</td>
      <td>4.850</td>
      <td>3.980295</td>
      <td>3.669143</td>
      <td>3.540811</td>
      <td>4.407405</td>
      <td>60.374169</td>
      <td>-24.634717</td>
      <td>7.328994</td>
      <td>1.601481</td>
      <td>5134765416778736512</td>
    </tr>
    <tr>
      <th>3</th>
      <td>57 Cet</td>
      <td>29.942590</td>
      <td>-20.824468</td>
      <td>5.429</td>
      <td>4.111464</td>
      <td>3.613912</td>
      <td>3.437497</td>
      <td>4.662973</td>
      <td>18.826217</td>
      <td>15.859864</td>
      <td>4.763460</td>
      <td>2.047709</td>
      <td>5136714576016018432</td>
    </tr>
    <tr>
      <th>4</th>
      <td>HD  12767</td>
      <td>31.122699</td>
      <td>-29.296752</td>
      <td>4.690</td>
      <td>5.099089</td>
      <td>4.642039</td>
      <td>4.677194</td>
      <td>4.641025</td>
      <td>1.849292</td>
      <td>11.772053</td>
      <td>9.403377</td>
      <td>-0.180284</td>
      <td>5021010046848175616</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>174</th>
      <td>47 Aqr</td>
      <td>335.398148</td>
      <td>-21.598607</td>
      <td>5.124</td>
      <td>4.504360</td>
      <td>4.247096</td>
      <td>4.141161</td>
      <td>4.810523</td>
      <td>-10.986031</td>
      <td>-84.833784</td>
      <td>17.347473</td>
      <td>1.235893</td>
      <td>6628587855877381120</td>
    </tr>
    <tr>
      <th>175</th>
      <td>bet PsA</td>
      <td>337.876686</td>
      <td>-32.346159</td>
      <td>4.290</td>
      <td>4.592107</td>
      <td>4.187389</td>
      <td>4.209581</td>
      <td>4.263408</td>
      <td>58.639817</td>
      <td>-18.857117</td>
      <td>22.083617</td>
      <td>0.027941</td>
      <td>6601750220152445440</td>
    </tr>
    <tr>
      <th>176</th>
      <td>rho Gru</td>
      <td>340.874958</td>
      <td>-41.414762</td>
      <td>4.835</td>
      <td>4.268914</td>
      <td>4.011698</td>
      <td>3.908653</td>
      <td>4.550108</td>
      <td>9.589843</td>
      <td>-93.616195</td>
      <td>13.415835</td>
      <td>1.187648</td>
      <td>6544701440869855744</td>
    </tr>
    <tr>
      <th>177</th>
      <td>HD 216149</td>
      <td>342.759098</td>
      <td>-39.156869</td>
      <td>5.400</td>
      <td>4.508897</td>
      <td>4.206003</td>
      <td>4.079573</td>
      <td>4.930226</td>
      <td>17.775482</td>
      <td>-7.583302</td>
      <td>3.627289</td>
      <td>1.580908</td>
      <td>6546508900546867840</td>
    </tr>
    <tr>
      <th>178</th>
      <td>HD 224554</td>
      <td>359.732844</td>
      <td>-52.745531</td>
      <td>5.133</td>
      <td>4.470931</td>
      <td>4.209739</td>
      <td>4.094971</td>
      <td>4.795882</td>
      <td>58.509812</td>
      <td>62.021589</td>
      <td>12.097462</td>
      <td>1.299198</td>
      <td>6521403820271516160</td>
    </tr>
  </tbody>
</table>
<p>179 rows × 13 columns</p>
</div>




```python
brighter_out = pd.DataFrame(columns=['simbad_name', 'ra', 'dec', 'Vmag', 'Imag','Sdss_mag_i', 'WFS6535 mag',
       'WFSHaIR mag', 'phot_g_mean_mag', 'pmra', 'pmdec', 'parallax',
       'bp_rp','source_id'])
brighter_out['simbad_name'] =  brighter['main_id']
brighter_out['ra'] = brighter['ra']
brighter_out['dec'] = brighter['dec']
brighter_out['Vmag'] = brighter['Vmag']
brighter_out['Imag'] = brighter['Imag']
brighter_out['pmra'] = brighter['pmra']
brighter_out['pmdec'] = brighter['pmdec']
for col in ['Sdss_mag_i', 'WFS6535 mag',
       'WFSHaIR mag', 'phot_g_mean_mag','parallax',
       'bp_rp','source_id']:
    brighter_out[col] = np.nan
```


```python
r_out_out = pd.concat([r_out, brighter_out], axis=0)
r_out_out = r_out_out.sort_values('ra')
r_out_out = r_out_out.reset_index(drop=True)
r_out_out
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>simbad_name</th>
      <th>ra</th>
      <th>dec</th>
      <th>Vmag</th>
      <th>Imag</th>
      <th>WFS6535 mag</th>
      <th>WFSHaIR mag</th>
      <th>phot_g_mean_mag</th>
      <th>pmra</th>
      <th>pmdec</th>
      <th>parallax</th>
      <th>bp_rp</th>
      <th>source_id</th>
      <th>Sdss_mag_i</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>30 Psc</td>
      <td>0.490081</td>
      <td>-6.014071</td>
      <td>NaN</td>
      <td>1.440000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>46.941000</td>
      <td>-40.471000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>1</th>
      <td>eps Phe</td>
      <td>2.352667</td>
      <td>-45.747424</td>
      <td>NaN</td>
      <td>2.590000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>120.393000</td>
      <td>-179.597000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>2</th>
      <td>7 Cet</td>
      <td>3.660070</td>
      <td>-18.932863</td>
      <td>NaN</td>
      <td>1.920000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-26.458000</td>
      <td>-73.450000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>3</th>
      <td>iot Cet</td>
      <td>4.856977</td>
      <td>-8.823919</td>
      <td>NaN</td>
      <td>2.110000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-14.610000</td>
      <td>-36.668000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>4</th>
      <td>alf Phe</td>
      <td>6.571048</td>
      <td>-42.305987</td>
      <td>2.380</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>233.050000</td>
      <td>-356.300000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>354</th>
      <td>phi Aqr</td>
      <td>348.580659</td>
      <td>-6.049006</td>
      <td>NaN</td>
      <td>1.860000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>36.575000</td>
      <td>-195.441000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>355</th>
      <td>psi01 Aqr</td>
      <td>348.972892</td>
      <td>-9.087735</td>
      <td>NaN</td>
      <td>2.900000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>369.477000</td>
      <td>-16.981000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>356</th>
      <td>b01 Aqr</td>
      <td>350.742608</td>
      <td>-20.100582</td>
      <td>NaN</td>
      <td>2.560000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-121.078000</td>
      <td>-97.781000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>357</th>
      <td>b02 Aqr</td>
      <td>351.511608</td>
      <td>-20.642018</td>
      <td>NaN</td>
      <td>2.460000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>-51.118000</td>
      <td>-63.995000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>358</th>
      <td>HD 224554</td>
      <td>359.732844</td>
      <td>-52.745531</td>
      <td>5.133</td>
      <td>4.470931</td>
      <td>4.209739</td>
      <td>4.094971</td>
      <td>4.795882</td>
      <td>58.509812</td>
      <td>62.021589</td>
      <td>12.097462</td>
      <td>1.299198</td>
      <td>6521403820271516160</td>
      <td>NaN</td>
    </tr>
  </tbody>
</table>
<p>359 rows × 14 columns</p>
</div>




```python
r_out_out['ra_hms'], r_out_out['dec_dms'] = np.nan, np.nan
from astropy.coordinates import SkyCoord
for i in range(len(r_out_out)):
    ob = SkyCoord(ra = r_out_out['ra'][i], dec = r_out_out['dec'][i], frame="icrs", unit="deg")
    # convert to string in hms and dms, and split the string in to [ra,dec]
    o = ob.to_string('hmsdms').split(' ')
    o = [o[i].replace('h',':') for i in [0,1]]
    o = [o[i].replace('m',':') for i in [0,1]]
    o = [o[i].replace('s','') for i in [0,1]]
    o = [o[i].replace('d',':') for i in [0,1]]
    # put into table:
    r_out_out.loc[i,'ra_hms'], r_out_out.loc[i,'dec_dms'] = o[0][:11],o[1][:12]

r_out_out['vetted'] = 'N'
r_out_out['n_sys'] = 1
r_out_out = r_out_out[['simbad_name','ra_hms', 'dec_dms','Vmag', 'Imag', 'WFS6535 mag',
       'WFSHaIR mag', 'phot_g_mean_mag', 'vetted', 'n_sys','pmra', 'pmdec', 'parallax',
       'bp_rp','source_id', 'ra', 'dec']]
r_out_out
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>simbad_name</th>
      <th>ra_hms</th>
      <th>dec_dms</th>
      <th>Vmag</th>
      <th>Imag</th>
      <th>WFS6535 mag</th>
      <th>WFSHaIR mag</th>
      <th>phot_g_mean_mag</th>
      <th>vetted</th>
      <th>n_sys</th>
      <th>pmra</th>
      <th>pmdec</th>
      <th>parallax</th>
      <th>bp_rp</th>
      <th>source_id</th>
      <th>ra</th>
      <th>dec</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>30 Psc</td>
      <td>00:01:57.61</td>
      <td>-06:00:50.65</td>
      <td>NaN</td>
      <td>1.440000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>N</td>
      <td>1</td>
      <td>46.941000</td>
      <td>-40.471000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.490081</td>
      <td>-6.014071</td>
    </tr>
    <tr>
      <th>1</th>
      <td>eps Phe</td>
      <td>00:09:24.64</td>
      <td>-45:44:50.72</td>
      <td>NaN</td>
      <td>2.590000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>N</td>
      <td>1</td>
      <td>120.393000</td>
      <td>-179.597000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>2.352667</td>
      <td>-45.747424</td>
    </tr>
    <tr>
      <th>2</th>
      <td>7 Cet</td>
      <td>00:14:38.41</td>
      <td>-18:55:58.30</td>
      <td>NaN</td>
      <td>1.920000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>N</td>
      <td>1</td>
      <td>-26.458000</td>
      <td>-73.450000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>3.660070</td>
      <td>-18.932863</td>
    </tr>
    <tr>
      <th>3</th>
      <td>iot Cet</td>
      <td>00:19:25.67</td>
      <td>-08:49:26.10</td>
      <td>NaN</td>
      <td>2.110000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>N</td>
      <td>1</td>
      <td>-14.610000</td>
      <td>-36.668000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>4.856977</td>
      <td>-8.823919</td>
    </tr>
    <tr>
      <th>4</th>
      <td>alf Phe</td>
      <td>00:26:17.05</td>
      <td>-42:18:21.55</td>
      <td>2.380</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>N</td>
      <td>1</td>
      <td>233.050000</td>
      <td>-356.300000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>6.571048</td>
      <td>-42.305987</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>354</th>
      <td>phi Aqr</td>
      <td>23:14:19.35</td>
      <td>-06:02:56.42</td>
      <td>NaN</td>
      <td>1.860000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>N</td>
      <td>1</td>
      <td>36.575000</td>
      <td>-195.441000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>348.580659</td>
      <td>-6.049006</td>
    </tr>
    <tr>
      <th>355</th>
      <td>psi01 Aqr</td>
      <td>23:15:53.49</td>
      <td>-09:05:15.84</td>
      <td>NaN</td>
      <td>2.900000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>N</td>
      <td>1</td>
      <td>369.477000</td>
      <td>-16.981000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>348.972892</td>
      <td>-9.087735</td>
    </tr>
    <tr>
      <th>356</th>
      <td>b01 Aqr</td>
      <td>23:22:58.22</td>
      <td>-20:06:02.09</td>
      <td>NaN</td>
      <td>2.560000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>N</td>
      <td>1</td>
      <td>-121.078000</td>
      <td>-97.781000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>350.742608</td>
      <td>-20.100582</td>
    </tr>
    <tr>
      <th>357</th>
      <td>b02 Aqr</td>
      <td>23:26:02.78</td>
      <td>-20:38:31.26</td>
      <td>NaN</td>
      <td>2.460000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>N</td>
      <td>1</td>
      <td>-51.118000</td>
      <td>-63.995000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>351.511608</td>
      <td>-20.642018</td>
    </tr>
    <tr>
      <th>358</th>
      <td>HD 224554</td>
      <td>23:58:55.88</td>
      <td>-52:44:43.91</td>
      <td>5.133</td>
      <td>4.470931</td>
      <td>4.209739</td>
      <td>4.094971</td>
      <td>4.795882</td>
      <td>N</td>
      <td>1</td>
      <td>58.509812</td>
      <td>62.021589</td>
      <td>12.097462</td>
      <td>1.299198</td>
      <td>6521403820271516160</td>
      <td>359.732844</td>
      <td>-52.745531</td>
    </tr>
  </tbody>
</table>
<p>359 rows × 17 columns</p>
</div>




```python
for i in range(len(r_out_out)):
    if np.isnan(r_out_out.loc[i,'pmra']):
        print(i)
    elif np.isnan(r_out_out.loc[i,'pmdec']):
        print(i)
```


```python
r_out_out = r_out_out.drop(index=[238,74,106,113,116,166,205,301,335])
r_out_out = r_out_out.reset_index(drop=True)
```


```python
r_out_out['Num'] = np.arange(1,len(r_out_out)+1,1)
r_out_out.to_csv('Bright-AO-Stars.csv', index=False)
```

## Plot of RA/Dec/RUWE of selected stars:


```python
%matplotlib inline
plt.scatter(r_out['ra'],r_out['dec'], c=r_out['phot_g_mean_mag'], cmap='viridis_r',vmin=min(r_out['phot_g_mean_mag']))
cbar = plt.colorbar()
cbar.ax.invert_yaxis()
cbar.ax.set_ylabel('Gaia g')
plt.xlabel('ra [deg]')
plt.ylabel('dec [deg]')
plt.grid(ls=':')
```


    
![png](images/output_60_0.png)
    


# Put catalog into TCS catalog format


```python
import warnings
warnings.filterwarnings('ignore')
pdcat = r_out_out.copy()
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
    pdcat['pmra s/yr'][i] = a3
    
    # Dec is easier:
    a = pdcat['pmdec'][i]*u.mas.to(u.arcsec)
    # put into table:
    pdcat['pmdec arcsec/yr'][i] = a
    
ind = np.argsort(pdcat['ra'])
pdcat = pdcat.loc[ind]
pdcat = pdcat.reset_index(drop=True)

pdcat['num'] = np.arange(1,len(pdcat)+1,1)
pdcat['Name'] = pdcat['simbad_name']

# pdcat['# Notes'] = np.nan
# for i in range(len(pdcat)):

#     pdcat['# Notes'][i] =\
#          '# WFS 65/35 mag: '+str(np.round(pdcat['WFS6535 mag'][i],decimals=2)) +\
#          ', WFS Ha/IR mag: '+str(np.round(pdcat['WFSHaIR mag'][i], decimals=2)) +\
#          ', Gaia g: '+str(np.round(pdcat['phot_g_mean_mag'][i], decimals=2)) +\
#          ', RUWE: '+str(pdcat['ruwe'][i])    

pdcat

```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>simbad_name</th>
      <th>ra_hms</th>
      <th>dec_dms</th>
      <th>Vmag</th>
      <th>Imag</th>
      <th>WFS6535 mag</th>
      <th>WFSHaIR mag</th>
      <th>phot_g_mean_mag</th>
      <th>vetted</th>
      <th>n_sys</th>
      <th>...</th>
      <th>parallax</th>
      <th>bp_rp</th>
      <th>source_id</th>
      <th>ra</th>
      <th>dec</th>
      <th>Num</th>
      <th>pmra s/yr</th>
      <th>pmdec arcsec/yr</th>
      <th>num</th>
      <th>Name</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>*  30 Psc</td>
      <td>00:01:57.61</td>
      <td>-06:00:50.65</td>
      <td>NaN</td>
      <td>1.440000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>N</td>
      <td>1</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>0.490081</td>
      <td>-6.014071</td>
      <td>1</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>1</td>
      <td>*  30 Psc</td>
    </tr>
    <tr>
      <th>1</th>
      <td>* eps Phe</td>
      <td>00:09:24.64</td>
      <td>-45:44:50.72</td>
      <td>NaN</td>
      <td>2.590000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>N</td>
      <td>1</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>2.352667</td>
      <td>-45.747424</td>
      <td>2</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>2</td>
      <td>* eps Phe</td>
    </tr>
    <tr>
      <th>2</th>
      <td>[KIS2018] A2744 46.3</td>
      <td>00:14:18.60</td>
      <td>-30:24:31.35</td>
      <td>NaN</td>
      <td>1.700000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>N</td>
      <td>1</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>3.577530</td>
      <td>-30.408710</td>
      <td>3</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>3</td>
      <td>[KIS2018] A2744 46.3</td>
    </tr>
    <tr>
      <th>3</th>
      <td>[KIS2018] A2744 46.1</td>
      <td>00:14:22.80</td>
      <td>-30:24:02.7</td>
      <td>NaN</td>
      <td>2.100000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>N</td>
      <td>1</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>3.595020</td>
      <td>-30.400750</td>
      <td>4</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>4</td>
      <td>[KIS2018] A2744 46.1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>*   7 Cet</td>
      <td>00:14:38.41</td>
      <td>-18:55:58.30</td>
      <td>NaN</td>
      <td>1.920000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>N</td>
      <td>1</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>3.660070</td>
      <td>-18.932863</td>
      <td>5</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>5</td>
      <td>*   7 Cet</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>355</th>
      <td>* phi Aqr</td>
      <td>23:14:19.35</td>
      <td>-06:02:56.42</td>
      <td>NaN</td>
      <td>1.860000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>N</td>
      <td>1</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>348.580659</td>
      <td>-6.049006</td>
      <td>356</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>356</td>
      <td>* phi Aqr</td>
    </tr>
    <tr>
      <th>356</th>
      <td>* psi01 Aqr</td>
      <td>23:15:53.49</td>
      <td>-09:05:15.84</td>
      <td>NaN</td>
      <td>2.900000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>N</td>
      <td>1</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>348.972892</td>
      <td>-9.087735</td>
      <td>357</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>357</td>
      <td>* psi01 Aqr</td>
    </tr>
    <tr>
      <th>357</th>
      <td>* b01 Aqr</td>
      <td>23:22:58.22</td>
      <td>-20:06:02.09</td>
      <td>NaN</td>
      <td>2.560000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>N</td>
      <td>1</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>350.742608</td>
      <td>-20.100582</td>
      <td>358</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>358</td>
      <td>* b01 Aqr</td>
    </tr>
    <tr>
      <th>358</th>
      <td>* b02 Aqr</td>
      <td>23:26:02.78</td>
      <td>-20:38:31.26</td>
      <td>NaN</td>
      <td>2.460000</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>N</td>
      <td>1</td>
      <td>...</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>351.511608</td>
      <td>-20.642018</td>
      <td>359</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>359</td>
      <td>* b02 Aqr</td>
    </tr>
    <tr>
      <th>359</th>
      <td>HD 224554</td>
      <td>23:58:55.88</td>
      <td>-52:44:43.91</td>
      <td>5.133</td>
      <td>4.470931</td>
      <td>4.209739</td>
      <td>4.094971</td>
      <td>4.795882</td>
      <td>N</td>
      <td>1</td>
      <td>...</td>
      <td>12.097462</td>
      <td>1.299198</td>
      <td>6521403820271516160</td>
      <td>359.732844</td>
      <td>-52.745531</td>
      <td>360</td>
      <td>0.003901</td>
      <td>0.062022</td>
      <td>360</td>
      <td>HD 224554</td>
    </tr>
  </tbody>
</table>
<p>360 rows × 22 columns</p>
</div>




```python
pdcat_out = pdcat[['num','Name','ra','dec']]
pdcat_out['Equinox'] = 2000.0
pdcat_out['pmra'] = pdcat['pmra s/yr']
pdcat_out['pmdec'] = pdcat['pmdec arcsec/yr'] 
pdcat_out['rotang'] = 0
pdcat_out['rot_mode'] = 'GRV'
pdcat_out['RA_probe1'],pdcat_out['Dec_probe1'] = '00:00:00.00',  '+00:00:00.0'
pdcat_out['equinox'] = 2000.0
pdcat_out['RA_probe2'],pdcat_out['Dec_probe2'] = '00:00:00.00',  '+00:00:00.0'
pdcat_out['equinox '] = 2000.0
pdcat_out['epoch'] = 2000.0
#pdcat_out['# Notes'] = pdcat['# Notes']

# pdcat_out.to_csv('Bright_AO_stars_cat.cat', index=False, sep='\t')
# pdcat_out.to_csv('Bright_AO_stars_cat.csv', index=False)
```


```python
pdcat_writeto = pdcat_out.to_csv(index=False, sep='\t')
pdcat_writeto
```

""")

