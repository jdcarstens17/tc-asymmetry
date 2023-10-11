# tc-asymmetry
This repository features MATLAB code to analyze tropical cyclone (TC) asymmetric structure (relative to vertical wind shear) in global reanalysis and climate model datasets. It includes generalized versions of scripts used to produce the figures in Carstens et al. (2023, J. Climate, submitted), as well as information about tropical cyclone tracks derived from the TempestExtremes software (https://github.com/ClimateGlobalChange/tempestextremes). It also includes scripts used to calculate vertical wind shear around each TC, and the corresponding wind shear data in .nc and .mat format. The scripts in this repository are designed to work with both .nc and .mat data, as described below. TC-centered snapshots of various kinematic and thermodynamic fields (as 20-degree-wide boxes) can be accessed through the Penn State Data Commons at <LINK PENDING>, and should be placed in the StormCenteredData folder for immediate use in a MATLAB workspace.

As of 11 October 2023: The repository accommodates 2012-2021 data from the ERA5 and CFSR reanalyses. StormCenteredData may neglect TC snapshots outside of the 0-30N latitude range by including NaN values outside this range, but this will be corrected within the next 4-6 weeks. Future releases will incorporate climate model datasets, including historical and warming scenario runs of the Community Atmosphere Model, into this same workflow, and new scripts will be added for in-depth process-oriented understanding of TC-wind shear interaction in these datasets.

# Description of Data - TCTracks
TCTracks contains files labeled "trajectories_<DATASET>.txt". These list track information for each tropical cyclone detected globally by TempestExtremes, ordered chronologically by track ID. Briefly, TempestExtremes searches for local minima in sea level pressure, co-located with a local maximum in geopotential thickness aloft to isolate warm core from cold core cyclones. It then stitches together snapshots in time to produce continuous tracks. Refer to Zarzycki and Ullrich (2017) for information on the particular parameters chosen to identify TCs using TempestExtremes.

Each line of the TC track file reads in the following order: track ID (from 0), year, month, day, hour, x and y index in the dataset corresponding to the TC center (from 0), longitude and latitude of that center, minimum sea level pressure (in Pa), maximum wind speed (in m/s). The data-processing scripts in this repository are built to load this information in immediately, then allow the user to choose whether to use wind speed or sea level pressure as a tool to place TC snapshots into composite bins. Other data files, including wind shear and the files in StormCenteredData, are structured in the exact same order as these track files. In other words, the 81st element in the time dimension for StormCenteredData corresponds to the 81st line of the track file. Scripts simply loop through all snapshots in order, using composite binning criteria set by the user to skip any snapshots as necessary.

# Description of Data - StormCenteredData
Data necessary to produce the figures and analysis in Carstens et al. (2023, J. Climate, submitted), as well as some supplemental fields, are included in this folder. The time dimension is simply the TC snapshot ID, in the same order as the TC track files described above. Otherwise, data are organized as 20-degree-wide boxes centered on the TC center identified by TempestExtremes (10 degrees in each direction). Data-processing scripts are already built to handle the differing horizontal resolutions of the datasets, where the user will only need to specify what dataset they would like to work with. No wind shear-relative rotation has been imposed at this stage. For vertically-varying fields such as wind and temperature, 27 vertical levels are captured from 100-1000 hPa. NaN is input for any missing data (i.e. data on the 1000 hPa pressure surface at the center of a 980 hPa TC).

While the data-processing scripts provided here are in MATLAB, the user has the option to work with either .mat or .nc-format data, depending on their preferred programming language. The naming convention for the files is <DATASET>_<VARIABLE>.nc, or <DATASET>_<VARIABLE>.mat. The available variables are described below:

|Variable Name|Description|
|---|---|
|`CIW`|Cloud ice water content (kg/kg). 4-D matrix (lon,lat,level,snapshot). NOT INCLUDED IN CFSR!|
|`CLW`|Cloud liquid water content in ERA5, combined cloud liquid and ice in CFSR (kg/kg). 4-D matrix (lon,lat,level,snapshot).|
|`CRR`|Precipitation rate (mm/hr) produced by the convective parameterization. 3-D matrix (lon,lat,snapshot).|
|`LHF`|Surface latent heat flux (W/m^2). 3-D matrix (lon,lat,snapshot).|
|`LSRR`|Precipitation rate (mm/hr) produced by the large-scale cloud scheme. 3-D matrix (lon,lat,snapshot).|
|`MLCAPE`|0-1 km mixed-layer convective available potential energy (CAPE, J/kg), calculated from T and Q fields. 3-D matrix (lon,lat,snapshot).|
|`MUCAPE`|Most-unstable parcel CAPE (J/kg), calculated from T and Q fields. 3-D matrix (lon,lat,snapshot).|
|`ModelCAPE`|Model-output CAPE (J/kg), which may differ between reanalyses! 3-D matrix (lon,lat,snapshot).|
|`Q`|Specific humidity (kg/kg). 4-D matrix (lon,lat,level,snapshot).|
|`SHF`|Surface sensible heat flux (W/m^2). 3-D matrix (lon,lat,snapshot).|
|`T`|Temperature (K). 4-D matrix (lon,lat,level,snapshot).|
|`U`|Zonal wind (m/s). 4-D matrix (lon,lat,level,snapshot).|
|`V`|Meridional wind (m/s). 4-D matrix (lon,lat,level,snapshot).|
|`W`|Vertical velocity in terms of pressure (Pa/s). Negative --> ascent! 4-D matrix (lon,lat,level,snapshot).|
|`Z`|Geopotential height (m). 4-D matrix (lon,lat,level,snapshot).|

# Description of Data - WindShear
Vertical wind shear information is included in a series of .mat and .nc files in this folder, ordered identically to the snapshots in StormCenteredData and TCTracks. The individual components of the shear are used to rotate snapshots in StormCenteredData in line with the shear, to allow assessment of shear-relative asymmetric TC structure. The file names contain information about the method of the shear calculation, as well as the vertical layer the calculation was performed over. The methods are: 1) Winds averaged over a 200-800 km radial annulus without vortex flow removal, 2) Winds averaged over a 0-500 km radial annulus without vortex flow removal, 3) Winds averaged over a 0-500 km radial annulus with the vortex flow removal technique of Davis et al. (2008, Mon. Wea. Rev.), and 4) Winds averaged over a 0-800 km radial annulus with the Davis et al. (2008) vortex flow removal. Standard deep-layer wind shear is defined as the vector difference between 850 and 200 hPa winds, while we also computed shear over shallower layers of 850-500 and 500-200 hPa.

A generalized script is included in this folder to calculate wind shear using the U and V winds in StormCenteredData, along with supplemental scripts for the vortex flow removal courtesy of Nicholas Barron. There, the user may select the vertical layer, radial coverage of the wind averaging, and whether or not to use the vortex removal. The resulting variables are described below as 1-D arrays:

|Variable Name|Description|
|---|---|
|`shear_u`|Zonal component of vertical wind shear (m/s).|
|`shear_v`|Meridional component of vertical wind shear (m/s).|
|`shear_magnitude`|Magnitude of vertical wind shear (m/s).|
|`shear_direction`|Azimuth of vertical wind shear (m/s).|

# Data-Processing Scripts
A brief description of each MATLAB script is below, taking in data from the StormCenteredData, TCTracks, and WindShear folders to produce wind shear-relative composite analyses of TC structure and processes. Each script is commented to describe the purpose, workflow, and resulting plots, with file paths aligned with the overall layout of this repository. Within the scripts, users have the ability to set their own binning schemes to assess composite TC structure, including ranges for TC wind speed, minimum pressure, and wind shear magnitude.

|Script|Purpose and Outputs|
|---|---|
|`CloudWater.m`|Displays cloud liquid and ice water content in radius-pressure space in 4 shear-relative quadrants. Also computes and overlays kinematic and thermodynamic properties on these plots, including vertical motion, convergence, and temperature fields.|
|`OmegaEquation.m`|Calculates and displays forcing terms for vertical motion using the generalized omega equation of Krishnamurti (1968, Mon. Wea. Rev.). This includes terms related to vorticity advection, vorticity tilting, and buoyancy advection, among others. Plots are shown in terms of spatial shear-relative maps of individual terms at each pressure level, and as radially-averaged vertical profiles, where users may set the radial boundaries.|
|`RainRate.m`|Plots shear-relative rainfall as both quadrant-specific radial profiles and spatial maps. The total rainfall is output, as well as contributions from the convective parameterization and large-scale cloud scheme. Other fields may also be overlaid, including vertical motion, convergence, and instability (CAPE). A section of code is also included (commented out) to perform a two-sided t-test for statistical significance between quadrant radial profiles of rainfall.|
|`SampleSizes.m`|Using binning schemes for TC intensity and wind shear set by the user, this script accumulates the sample sizes of each composite bin, and outputs them as bar graphs.|
|`TanRadVertWinds.m`|Calculates and plots tangential and radial winds at different vertical levels specified by the user in shear-relative space, overlaying the vertical motion field as well. Spatial maps and radial profiles are produced similarly to RainRate.m, where code is included to perform two-sided t-test statistical significance testing between quadrant radial profiles of tangential and radial wind.|
|`Thermo.m`|Calculates and displays various thermodynamic fields in shear-relative space, both in terms of raw values and as anomalies from the azimuthal mean. This includes specific humidity, relative humidity, potential temperature, and equivalent potential temperature, where users may overlay vertical motion fields and specify the vertical level of interest, as well as surface latent and sensible heat fluxes.|
|`VortexTilt.m`|Calculates and plots composite mean vortex tilt, where the TC center is identified above the surface in 100 hPa increments as either a centroid of cyclonic relative vorticity (default), or a centroid of negative geopotential height anomalies. This is done from 900-400 hPa.|
|`flowfun.m`|Utility script to calculate velocity potential and streamfunction, called within OmegaEquation.m (originally written by Kirill Pankratov).|
|`simpson_summation.m`|Utility script to perform Simpson-rule column-wise cumulative summation, called within flowfun.m (and by extension, OmegaEquation.m) (originally written by Kirill Pankratov).|
|`tanradwinds.m`|Utility script to calculate tangential and radial winds from the zonal and meridional components, for one vertical level at a time.|
|`lldistkm.m`|Utility script to calculate the distance between two points in km, from lat/lon coordinates. Available for download at https://www.mathworks.com/matlabcentral/fileexchange/38812-latlon-distance. Copyright (c) 2012, M Sohrabinia, All rights reserved. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.|
