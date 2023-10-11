% FUNCTION Calc_WindShear - Edited by Jake Carstens 10/11/2023 for manuscript on tropical cyclone asymmetry in reanalyses
%                           Generalized version of script used to generate wind shear data files for all analysis.
% PURPOSE: Compute vertical wind shear over a vertical layer and radius range specified by the user, with the option
%          to employ the vortex removal strategy of Davis et al. (2008), subtracting irrotational and nondivergent wind.
% NOTES: 1. Storm-centered data snapshots and TempestExtremes-derived TC tracks will be provided in same directory
%           within repository, so filepaths should not need to be altered.
%        2. The user has the choice to use either .mat or .nc format for the storm-centered data. TempestExtremes
%           tracks are provided as .txt files.
% PROCEDURE:
% 1. Load in TC tracks and information about the grid.
% 2. Pre-allocate matrices for u/v components of wind shear, as well as magnitude and direction.
% 3. Loop through all TC snapshots, select wind shear calculation strategy, and perform calculation. Save as both .mat and .nc files.

% FIRST, STATE WHAT REANALYSIS/MODEL GRID WE'RE USING TO DEFINE THE 10-DEGREE BOXES
dset="CFSR";
if (dset == "CFSR")
  x10d=10; y10d=10; xtotal=720; xgrid=0.50; ygrid=0.50; dset_c='CFSR';
elseif (dset == "ERA5")
  x10d=20; y10d=20; xtotal=1440; xgrid=0.25; ygrid=0.25; dset_c='ERA5';
end
lev=[100:25:250 300:50:750 775:25:1000];   % 27 total vertical levels.

% THEN, LIST THE DESIRED VERTICAL LEVELS AND RADII THAT YOU WISH TO PERFORM THE CALCULATION OVER.
% LIST VERTICAL LEVELS AS INDICIES OF lev --> "7" refers to 250 hPa, for example.
% LIST RADII IN KM. CURRENTLY SET FOR 850-200 hPa LAYER OVER A 0-500 KM ANNULUS.
% YOU WILL ALSO BE ABLE TO CHOOSE WHETHER OR NOT TO EMPLOY VORTEX REMOVAL HERE.
lev_lower=21; lev_upper=5; innerbound=0; outerbound=500;
vortex_removal="ON"; <-- currently set to remove vortex flow.
if (vortex_removal == "ON")   % This will help automate the resulting file name to reflect how the shear was calculated.
  method_text='VortexRemoval';
else
  method_text='Annulus';
end

% LOAD IN TC TRACKS FROM TEMPESTEXTREMES. LOAD IN LAT/LON TO IDENTIFY HOW LONG WE NEED OUR WIND SHEAR ARRAYS TO BE.
A=readmatrix(['TCTracks/trajectories_' dset_c '.txt']); clon=A(:,8); clat=A(:,9); clear A

% INITIALIZE ARRAYS FOR U/V SHEAR COMPONENTS, MAGNITUDE, AND DIRECTION.
shear_u=[]; shear_v=[]; shear_magnitude=[]; shear_direction=[];

% IF WORKING WITH .mat DATA, QUEUE UP THE FILES FOR ACCESS HERE. IF CHOOSING .nc, COMMENT THESE LINES OUT.
U=matfile(['StormCenteredData/' dset_c '_U.mat']); V=matfile(['StormCenteredData/' dset_c '_V.mat']);

% ENTER LOOP MOVING THROUGH ALL TC SNAPSHOTS, WHERE RELEVANT DATA WILL BE LOADED, ROTATED ABOUT THE SHEAR, AND PLACED INTO COMPOSITES.
for p=1:length(clon)
  tic                              % Matlab utility, along with "toc", to track how long it takes to run each iteration of the loop.
% LOAD DATA HERE! BE SURE TO COMMENT OUT WHATEVER DATA FORMAT YOU ARE NOT USING BETWEEN .nc AND .mat!
% .mat SECTION
  u=U.u(:,:,:,p); v=V.v(:,:,:,p);

% .nc SECTION
%  start=[1 1 1 p]; interval=[length([-10:xgrid:10]) length([-10:ygrid:10]) length(lev) 1];
%  u=ncread(['StormCenteredData/' dset_c '_U.nc'],'u',start,interval);
%  v=ncread(['StormCenteredData/' dset_c '_V.nc'],'v',start,interval);

  U_lower=u(:,:,lev_lower);
  V_lower=v(:,:,lev_lower);
  U_upper=u(:,:,lev_upper);    % Set lev_lower and lev_upper above, when first loading in the TC tracks.
  V_upper=v(:,:,lev_upper);

  dx=lldistkm([clat(p) clon(p)],[clat(p) clon(p)+xgrid]).*1000.0;   % Approximate zonal grid spacing, converted to m.
  dy=lldistkm([clat(p) clon(p)],[clat(p)+ygrid clon(p)]).*1000.0;   % Approximate meridional grid spacing, converted to m.
  lat=[-10:ygrid:10]+clat(p); lon=[-10:xgrid:10]+clon(p);           % Set storm-relative latitude/longitude arrays

  dist=[];   % Create matrix of distances between each point in box and TC center
  for x=1:length(lon)
    for y=1:length(lat)
      dist(x,y)=lldistkm([clat(p) clon(p)],[lat(y) lon(x)]);
    end
  end
  [uselon,uselat]=find(dist>=innerbound & dist<=outerbound);  % Return indices where our annulus radial condition is met.

  if (vortex_removal == "ON")
    % Calculate nondivergent and irrotational winds throughout the 20-degree TC-centered box at both vertical levels.
    % NEED TO HAVE THE FOLLOWING UTILITIES IN WORKING DIRECTORY: helmholtz_decompose.m, d_dx_c.m
    [u_psi_lower,v_psi_lower,u_chi_lower,v_chi_lower,psi_lower,chi_lower]=helmholtz_decompose(U_lower,V_lower,dx,dy);
    [u_psi_upper,v_psi_upper,u_chi_upper,v_chi_upper,psi_upper,chi_upper]=helmholtz_decompose(U_upper,V_upper,dx,dy);

    % Remove the nondivergent and irrotational winds everywhere but the boundary.
    U_lower(2:end-1,2:end-1) = U_lower(2:end-1, 2:end-1)-u_psi_lower(2:end-1,2:end-1)-u_chi_lower(2:end-1,2:end-1);
    V_lower(2:end-1,2:end-1) = V_lower(2:end-1, 2:end-1)-v_psi_lower(2:end-1,2:end-1)-v_chi_lower(2:end-1,2:end-1);
    U_upper(2:end-1,2:end-1) = U_upper(2:end-1, 2:end-1)-u_psi_upper(2:end-1,2:end-1)-u_chi_upper(2:end-1,2:end-1);
    V_upper(2:end-1,2:end-1) = V_upper(2:end-1, 2:end-1)-v_psi_upper(2:end-1,2:end-1)-v_chi_upper(2:end-1,2:end-1);
  end

  % ACCUMULATE RESULTING U/V WINDS INTO 1-D ARRAYS TO USE IN CALCULATING MEAN FLOW AND SHEAR.
  useU_lower=[]; useU_upper=[]; useV_lower=[]; useV_upper=[];
  for a=1:length(uselat)
    useU_lower(a)=U_lower(uselon(a),uselat(a));
    useU_upper(a)=U_upper(uselon(a),uselat(a));
    useV_lower(a)=V_lower(uselon(a),uselat(a));
    useV_upper(a)=V_upper(uselon(a),uselat(a));
  end

  % CALCULATE U AND V COMPONENTS OF WIND SHEAR, AND ADD TO EXISTING U, V, MAGNITUDE, AND DIRECTION ARRAYS
  shear_u(p)=mean(useU_upper)-mean(useU_lower);
  shear_v(p)=mean(useV_upper)-mean(useV_lower);
  shear_magnitude(p)=sqrt(shear_u(p)^2+shear_v(p)^2); %in m/s
  sheardir=atand(shear_v(p)/shear_u(p)); %does arctan in degrees of v/u

  % REWORK SHEAR DIRECTION TO REFLECT DEGREES CLOCKWISE FROM POINTING DUE NORTH.
  if (shear_u(p) > 0 && shear_v(p) > 0)     % Upper right quadrant
    shear_direction(p) = 270-sheardir;
  elseif (shear_u(p) < 0 && shear_v(p) > 0) % Upper left quadrant
    shear_direction(p) = 90-sheardir;
  elseif (shear_u(p) < 0 && shear_v(p) < 0) % Lower left quadrant
    shear_direction(p) = 90-sheardir;
  else                                      % Lower right quadrant
    shear_direction(p) = 270-sheardir;
  end
  toc
end

% SAVE OUR 4 VARIABLES TO BOTH A .MAT FILE FOR EASY MATLAB ACCESSIBILITY...
matfile=['WindShear/Shear_' method_text '_' num2str(innerbound) 'to' num2str(outerbound) 'km_' dset_c '.mat'];
save(matfile,'shear_magnitude','shear_direction','shear_u','shear_v');

% ...AND TO A NETCDF FILE FOR USE WITH PYTHON OR OTHER LANGUAGES.
ncfile=['WindShear/Shear_' method_text '_' num2str(innerbound) 'to' num2str(outerbound) 'km_' dset_c '.nc'];
nccreate(ncfile,'snapshot','Datatype','int32','Dimensions',{'snapshot',length(clat)});
ncwrite(ncfile,'snapshot',1:length(clat));
nccreate(ncfile,'shear_u','Datatype','single','Dimensions',{'snapshot',length(clat)});
ncwrite(ncfile,'shear_u',shear_u);
nccreate(ncfile,'shear_v','Datatype','single','Dimensions',{'snapshot',length(clat)});
ncwrite(ncfile,'shear_v',shear_v);
nccreate(ncfile,'shear_magnitude','Datatype','single','Dimensions',{'snapshot',length(clat)});
ncwrite(ncfile,'shear_magnitude',shear_magnitude);
nccreate(ncfile,'shear_direction','Datatype','single','Dimensions',{'snapshot',length(clat)});
ncwrite(ncfile,'shear_direction',shear_direction);

ncwriteatt(ncfile,'snapshot','long_name','Snapshot ID');
ncwriteatt(ncfile,'shear_u','long_name','Zonal Wind Shear');
ncwriteatt(ncfile,'shear_v','long_name','Meridional Wind Shear');
ncwriteatt(ncfile,'shear_magnitude','long_name','Wind Shear Magnitude');
ncwriteatt(ncfile,'shear_direction','long_name','Wind Shear Direction');
ncwriteatt(ncfile,'snapshot','units',' ');
ncwriteatt(ncfile,'shear_u','units','m/s');
ncwriteatt(ncfile,'shear_v','units','m/s');
ncwriteatt(ncfile,'shear_magnitude','units','m/s');
ncwriteatt(ncfile,'shear_direction','units','Degrees Clockwise from Due North');
