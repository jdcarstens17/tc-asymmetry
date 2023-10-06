% FUNCTION OmegaEquation - Edited by Jake Carstens 10/6/2023 for manuscript on tropical cyclone asymmetry in reanalyses
%                          Generalized version of script used to generate Figures 9-10 of manuscript.
% PURPOSE: Calculate terms of generalized omega equation (Krishnamurti 1968, Zhang et al. 2000, Yu and Didlake 2019)
%          in shear-relative quadrants, and plot vertically-resolved heatmaps and single-level spatial maps.
% NOTES: 1. Storm-centered data snapshots, TempestExtremes-derived TC tracks, and wind shear files will be provided
%           in same directory within repository, so filepaths should not need to be altered.
%        2. The user has the choice to use either .mat or .nc format for the wind shear and storm-centered data. TempestExtremes
%           tracks are provided as .txt files.
%        3. REQUIRES EXTERNAL UTILITY FUNCTIONS flowfun.m AND simpson_summation.m (Copyright Kirill Pankratov)
% PROCEDURE:
% 1. Choose reanalysis dataset (currently ERA5 or CFSR), loading in its TC tracks, associated wind shear,
%    and information about the grid. Here, select vertical levels for vertically-resolved analysis.
% 2. Pre-allocate matrices for Laplacian of omega (v1), differential vertical advection of vorticity (v2),
%    differential vortex tilting (v3), differential horizontal advection of vorticity (v4), buoyancy advection (v5),
%    and a residual (mainly diabatic) term (v6), binned by both TC behavior and wind shear magnitude.
% 3. Loop through all TC snapshots, rotate data to shear vector, and place into TC/shear bin composite.
% 4. For each composite bin, plot single-level spatial maps of each term, as well as quadrant-specific heatmaps
%    encompassing all vertical levels chosen for analysis at the start. Positive values (red) will be considered as
%    contributing factors to ascent. Choose radial range to average heatmaps over within the plotting section.

% FIRST, CHOOSE THE DATASET. THIS INCLUDES INFORMATION ABOUT THE GRID SPACING, USEFUL FOR DRAWING A 10-DEGREE BOX
% AROUND THE TC (x10d/y10d), AS WELL AS RE-CENTERING (xtotal) THE DOMAIN ABOUT THE TC.
dset="CFSR";
if (dset == "CFSR")
  x10d=10; y10d=10; xtotal=720; xgrid=0.50; ygrid=0.50; dset_c='CFSR';
elseif (dset == "ERA5")
  x10d=20; y10d=20; xtotal=1440; xgrid=0.25; ygrid=0.25; dset_c='ERA5';
end
% Current structure is to analyze between 950-300 hPa, with 975 and 250 as bounds
plevs=[100:25:250 300:50:750 775:25:1000]; levs=8:1:25; plotlevs=plevs(levs);

% LOAD IN TC TRACKS FROM TEMPESTEXTREMES. LOAD IN STORM ID, TIME, POSITION, WIND SPEED, AND MINIMUM PRESSURE.
A=readmatrix(['TCTracks/trajectories_' dset_c '.txt']);  % Change to appropriate directory
id=single(A(:,1)+1); year=single(A(:,2)); month=single(A(:,3)); day=single(A(:,4)); hour=single(A(:,5));
i=single(A(:,6)+1); j=single(A(:,7)+1); lon=A(:,8); lat=A(:,9); pres=A(:,10); wind=A(:,11);
clear A

% LOAD IN WIND SHEAR DATA, AND SET BIN VALUES THAT WILL BE USED FOR COMPOSITING AND PLOT TITLES.
% OTHER WIND SHEAR FILES WILL BE INCLUDED, INCLUDING RESULTS USING DIFFERENT VERTICAL LAYERS AND A VORTEX REMOVAL STRATEGY.
load(['WindShear/Shear_Annulus_200to800km_' dset_c '.mat']); shearvar=shear_magnitude;    % .mat OPTION
%shearvar=ncread(['WindShear/Shear_Annulus_200to800km_' dset_c '.nc'],'shear_magnitude');  % .nc OPTION
%shear_u=ncread(['WindShear/Shear_Annulus_200to800km_' dset_c '.nc'],'shear_u');
%shear_v=ncread(['WindShear/Shear_Annulus_200to800km_' dset_c '.nc'],'shear_v');
shearbins=[0:5:10 30];     % Bin values for compositing. Shear is in m/s.
shearbintext=["Low Shear (< 5 m/s)","Mod. Shear (5-10 m/s)","High Shear (> 10 m/s)"];  % This will go in the plot titles.
shearplotname="ShrMag";    % This will go in the plot file names.

% DECIDE ON A SCHEME TO BIN THE TCs. OPTIONS INCLUDE WIND SPEED, MINIMUM PRESSURE, AND WIND SPEED CHANGE IN THE NEXT 6 HOURS.
tcbinscheme=2;
if (tcbinscheme == 1)
  tcbins=[10 17 25 33 100]; tcbintext=["TD","Weak TS","Strong TS","Hurricane"];  % Winds are in m/s. Major hurricanes are quite rare.
  tcvar=wind;          % Extracts maximum wind speed from each TempestExtremes snapshot.
  tcplotname="VMax";   % This will go in the plot file names, like shearplottitle above.
elseif (tcbinscheme == 2)
  tcbins=[900 980:10:1010]; tcbintext=["< 980 mb","980-990 mb","990-1000 mb","> 1000 mb"];
  tcvar=pres./100.0;   % Pressure is in Pa in TempestExtremes, so convert to hPa here.
  tcplotname="PMin";
elseif (tcbinscheme == 3)
  tcbins=[-10 -2 2 10];  % 2 m/s change over 6 hours equates to about 10 knots in 24 h.
  tcbintext=["Weakening","Steady-State","Intensifying"];
  tcvar=[]; tcplotname="Rate";  % Leave the TC variable blank for now, because I'll need to check each time to ensure the track doesn't end.
end

% INITIALIZE COMPOSITE MATRICES, INCLUDING A COUNT MATRIX FOR THE SNAPSHOTS FALLING INTO EACH TC, SHEAR, AND RADIUS BIN.
% TERMS ARE AS FOLLOWS, AND WILL BE COMMENTED IN EACH OF THEIR RESPECTIVE CALCULATION SECTIONS:
% 1. Laplacian of omega
% 2. Differential vertical advection of vorticity
% 3. Differential vortex tilting
% 4. Differential horizontal advection of vorticity
% 5. Buoyancy advection
% 6. Residual term, assumed to be primarily due to diabatic effects
dist=[0:25:500];                                                        % Radii for quadrant profiles will be rounded to nearest 25 km (~0.25 deg).
count=zeros(length(levs),length(tcbins)-1,length(shearbins)-1);                      % Accumulate all snapshots in a particular TC/shear bin.
v1_comp=zeros(2*x10d+1,2*y10d+1,length(levs),length(tcbins)-1,length(shearbins)-1);  % Draw a 10-deg by 10-deg box for the TC composite.
v2_comp=zeros(2*x10d+1,2*y10d+1,length(levs),length(tcbins)-1,length(shearbins)-1);  % Do this for each of "levs" specified above.
v3_comp=zeros(2*x10d+1,2*y10d+1,length(levs),length(tcbins)-1,length(shearbins)-1);  % Make sure that "levs" does not include the endpoints of the
v4_comp=zeros(2*x10d+1,2*y10d+1,length(levs),length(tcbins)-1,length(shearbins)-1);  % vertical level array! We need these for centered finite
v5_comp=zeros(2*x10d+1,2*y10d+1,length(levs),length(tcbins)-1,length(shearbins)-1);  % differencing bounds.
v6_comp=zeros(2*x10d+1,2*y10d+1,length(levs),length(tcbins)-1,length(shearbins)-1);
v1_quad=NaN(4,length(dist),length(levs),length(tcbins)-1,length(shearbins)-1,2500);  % Quadrant radial profiles. Order: DR --> DL --> UL --> UR
v2_quad=NaN(4,length(dist),length(levs),length(tcbins)-1,length(shearbins)-1,2500);  % I will keep these discrete and use "nanmean" to calculate
v3_quad=NaN(4,length(dist),length(levs),length(tcbins)-1,length(shearbins)-1,2500);  % the quadrant averages later.
v4_quad=NaN(4,length(dist),length(levs),length(tcbins)-1,length(shearbins)-1,2500);  % For other analyses, I have used this approach to allow for
v5_quad=NaN(4,length(dist),length(levs),length(tcbins)-1,length(shearbins)-1,2500);  % statistical significance testing between quadrant pairs.
v6_quad=NaN(4,length(dist),length(levs),length(tcbins)-1,length(shearbins)-1,2500);

% IF CHOOSING THE .mat OPTION, QUEUE UP THE FILES HERE
U=matfile(['StormCenteredData/' dset_c '_U.mat']);
V=matfile(['StormCenteredData/' dset_c '_V.mat']);
W=matfile(['StormCenteredData/' dset_c '_W.mat']);
T=matfile(['StormCenteredData/' dset_c '_T.mat']);

% LOOP THROUGH ALL TC SNAPSHOTS, WHERE TERMS WILL BE CALCULATED, ROTATED ABOUT THE SHEAR, AND PLACED INTO COMPOSITES.
for p=1:length(year)
  if (lat(p) > 30 || lat(p) < 0)   % Exclude extratropical and Southern Hemisphere snapshots.
    continue;                      % If SH snapshots are included, data will need to be mirrored across shear vector!
  end
  if (tcbinscheme == 3)            % Compute 6-hour intensification rate here, or skip ahead if it's a TC's last timestep
    if (id(p+1) ~= id(p))
      continue;
    elseif (p == length(year))     % Inherently, if we've reached the last timestep we should break out of this loop.
      break;
    else
      tcvar(p)=wind(p+1)-wind(p);
    end
  end

  tic                              % Matlab utility, along with "toc", to track how long it takes to run each iteration of the loop.

% LOAD DATA HERE! CHANGE FILE DIRECTORY AND NAMING STRUCTURE ACCORDINGLY IF REPRODUCING.
% .mat SECTION
  v_1=U.u(:,:,:,p); % U wind
  v_2=V.v(:,:,:,p); % V wind
  v_3=W.w(:,:,:,p); % W wind
  v_5=T.temp(:,:,:,p); % Temperature

% .nc SECTION
  start=[1 1 1 p]; interval=[length([-10:xgrid:10]) length([-10:ygrid:10]) length(lev) 1];
  v_1=ncread(['StormCenteredData/' dset_c '_U.nc'],'u',start,interval);
  v_2=ncread(['StormCenteredData/' dset_c '_V.nc'],'v',start,interval);
  v_3=ncread(['StormCenteredData/' dset_c '_W.nc'],'w',start,interval);
  v_5=ncread(['StormCenteredData/' dset_c '_T.nc'],'temp',start,interval);

  v_4=zeros(4*x10d+1,4*y10d+1,length(plevs)); % Relative vorticity
  for w=1:length(plevs)
    for x=2:4*x10d
      for y=2:4*y10d
        dx=lldistkm([latbox(y) lonbox(x-1)],[latbox(y) lonbox(x+1)]).*1000.0;
        dy=lldistkm([latbox(y-1) lonbox(x)],[latbox(y+1) lonbox(x)]).*1000.0;
        du=v_1(x,y+1,w)-v_1(x,y-1,w);
        dv=v_2(x+1,y,w)-v_2(x-1,y,w);
        v_4(x,y,w)=(dv./dx)-(du./dy);
      end
    end
  end

% CALCULATE TERMS OF OMEGA EQUATION HERE. REFER TO EQUATION IN YU AND DIDLAKE (2019) OR CARSTENS ET AL. (2023) DRAFT
% FOR THE TERMS INCLUDED HERE. I HAVE ATTEMPTED TO LAY OUT EACH CALCULATION STEP BY STEP USING VARIABLE NAMES MATCHING
% THE SYMBOLS IN THE EQUATION. NOTE THAT I HAVE LOADED IN A LARGER TC-CENTERED BOX OF DATA THAN WHAT WILL BE PLOTTED,
% AND AVOIDED USING THE VERTICAL BOUNDARIES IN MY "levs" ANALYSIS ARRAY. THIS IS BECAUSE MULTIPLE CENTERED FINITE
% DIFFERENCES WILL BE CALCULATED, REQUIRING LARGER 3-D BOUNDARIES. THIS SECTION REQUIRES THE FOLLOWING DOWNLOADED
% UTILITIES IN THE SAME FOLDER AS THIS FUNCTION: lldistkm, flowfun, simpson_summation

% -----------------------------------------------------------------------
% V1: LAPLACIAN OF OMEGA

  for z=1:length(plevs)
    theta(:,:,z)=v_5(:,:,z).*((1000./plevs(z)).^0.286);
    i_rho(:,:,z)=287.*v_5(:,:,z)./(plevs(z).*100);
  end
  for z=1:length(levs)
    dtdp(:,:,z)=(theta(:,:,levs(z)+1)-theta(:,:,levs(z)-1))./((plevs(levs(z)+1)-plevs(levs(z)-1)).*100.0);
  end
  sigma=-i_rho(:,:,levs)./theta(:,:,levs).*dtdp;
  sigma_w=sigma.*v_3(:,:,levs);
  cor=zeros(4*x10d+1,4*y10d+1,length(levs)); lat_box=lat(p)+latbox;
  for a=1:length(lat_box)
    cor(:,a,:)=2.*7.292e-5.*sind(lat_box(a));
  end
  for z=2:length(plevs)-1
    dwdp(:,:,z)=(v_3(:,:,z+1)-v_3(:,:,z-1))./((plevs(z+1)-plevs(z-1)).*100.0);
  end
  for z=levs
    dwdp2(:,:,z-min(levs)+1)=(dwdp(:,:,z+1)-dwdp(:,:,z-1))./((plevs(z+1)-plevs(z-1)).*100.0);
  end
  grad_w_x=zeros(4*x10d+1,4*y10d+1,length(levs));
  grad_w_y=zeros(4*x10d+1,4*y10d+1,length(levs));
  grad2_w=zeros(4*x10d+1,4*y10d+1,length(levs));
  for i=2:4*x10d
    for j=2:4*y10d
      dx=lldistkm([lat(p)+latbox(j) lon(p)+lonbox(i-1)],[lat(p)+latbox(j) lon(p)+lonbox(i+1)]).*1000.0; dy=lldistkm([lat(p)+latbox(j-1) lon(p)+lonbox(i)],$
      grad_w_x(i,j,:)=(sigma_w(i+1,j,:) - sigma_w(i-1,j,:))./dx;
      grad_w_y(i,j,:)=(sigma_w(i,j+1,:) - sigma_w(i,j-1,:))./dy;
    end
  end
  for i=3:4*x10d-1
    for j=3:4*y10d-1
      dx=lldistkm([lat(p)+latbox(j) lon(p)+lonbox(i-1)],[lat(p)+latbox(j) lon(p)+lonbox(i+1)]).*1000.0; dy=lldistkm([lat(p)+latbox(j-1) lon(p)+lonbox(i)],$
      grad2_w(i,j,:)=(grad_w_x(i+1,j,:) - grad_w_x(i-1,j,:))./dx + (grad_w_y(i,j+1,:) - grad_w_y(i,j-1,:))./dy;
    end
  end
  v1=zeros(4*x10d+1,4*y10d+1,length(levs));
  v1=grad2_w+(cor.^2.*dwdp2);
  clear dtdp sigma sigma_w dwdp dwdp2 grad_w_x grad_w_y grad2_w

% -----------------------------------------------------------------------
% V2: DIFFERENTIAL VERTICAL ADVECTION OF VORTICITY

  for z=1:length(plevs)
    psi(:,:,z)=flowfun(v_1(:,:,z),v_2(:,:,z),'-');
  end
  grad_psix=zeros(4*x10d+1,4*y10d+1,length(plevs));
  grad_psiy=zeros(4*x10d+1,4*y10d+1,length(plevs));
  grad2_psi=zeros(4*x10d+1,4*y10d+1,length(plevs));
  gradv_psi=zeros(4*x10d+1,4*y10d+1,length(plevs));
  for i=2:4*x10d
    for j=2:4*y10d
      dx=lldistkm([lat(p)+latbox(j) lon(p)+lonbox(i-1)],[lat(p)+latbox(j) lon(p)+lonbox(i+1)]).*1000.0; dy=lldistkm([lat(p)+latbox(j-1) lon(p)+lonbox(i)],$
      grad_psix(i,j,:)=(psi(i+1,j,:) - psi(i-1,j,:))./dx;
      grad_psiy(i,j,:)=(psi(i,j+1,:) - psi(i,j-1,:))./dy;
    end
  end
  for i=3:4*x10d-1
    for j=3:4*y10d-1
      dx=lldistkm([lat(p)+latbox(j) lon(p)+lonbox(i-1)],[lat(p)+latbox(j) lon(p)+lonbox(i+1)]).*1000.0; dy=lldistkm([lat(p)+latbox(j-1) lon(p)+lonbox(i)],$
      grad2_psi(i,j,:)=(grad_psix(i+1,j,:) - grad_psix(i-1,j,:))./dx + (grad_psiy(i,j+1,:) - grad_psiy(i,j-1,:))./dy;
    end
  end
  for z=2:length(plevs)-1
    gradv_psi(:,:,z)=(grad2_psi(:,:,z+1)-grad2_psi(:,:,z-1))./((plevs(z+1)-plevs(z-1)).*100.0).*v_3(:,:,z);
  end
  v2=zeros(4*x10d+1,4*y10d+1,length(levs));
  for z=levs
    v2(:,:,z-min(levs)+1)=(gradv_psi(:,:,z+1)-gradv_psi(:,:,z-1))./((plevs(z+1)-plevs(z-1)).*100.0).*cor(:,:,z-min(levs)+1);
  end
  clear grad_psix grad_psiy grad2_psi gradv_psi

% -----------------------------------------------------------------------
% V3: DIFFERENTIAL VORTEX TILTING

  dpsidp=zeros(4*x10d+1,4*y10d+1,length(plevs));
  dpsidx=zeros(4*x10d+1,4*y10d+1,length(plevs));
  dpsidy=zeros(4*x10d+1,4*y10d+1,length(plevs));
  dwdx=zeros(4*x10d+1,4*y10d+1,length(plevs)); dwdy=zeros(4*x10d+1,4*y10d+1,length(plevs));
  for z=2:length(plevs)-1
    dpsidp(:,:,z)=(psi(:,:,z+1)-psi(:,:,z-1))./((plevs(z+1)-plevs(z-1)).*100.0);
  end
  for i=2:4*x10d
    for j=2:4*y10d
      dx=lldistkm([lat(p)+latbox(j) lon(p)+lonbox(i-1)],[lat(p)+latbox(j) lon(p)+lonbox(i+1)]).*1000.0; dy=lldistkm([lat(p)+latbox(j-1) lon(p)+lonbox(i)],$
      dwdx(i,j,:)=(v_3(i+1,j,:) - v_3(i-1,j,:))./dx; dwdy(i,j,:)=(v_3(i,j+1,:) - v_3(i,j-1,:))./dy;
      dpsidx(i,j,:)=(dpsidp(i+1,j,:) - dpsidp(i-1,j,:))./dx; dpsidy(i,j,:)=(dpsidp(i,j+1,:) - dpsidp(i,j-1,:))./dy;
    end
  end
  tilting=(dwdx.*dpsidx)+(dwdy.*dpsidy);
  v3=zeros(4*x10d+1,4*y10d+1,length(levs));
  for z=levs
    v3(:,:,z-min(levs)+1)=(tilting(:,:,z+1)-tilting(:,:,z-1))./((plevs(z+1)-plevs(z-1)).*100.0).*cor(:,:,z-min(levs)+1);
  end
  clear dpsidp dpsidx dpsidy dwdx dwdy psi tilting

% -----------------------------------------------------------------------
% V4: DIFFERENTIAL VORTICITY ADVECTION

  vor_adv_u=zeros(4*x10d+1,4*y10d+1,length(plevs)); vor_adv_v=zeros(4*x10d+1,4*y10d+1,length(plevs)); vor_adv=zeros(4*x10d+1,4*y10d+1,length(plevs));
  for i=2:4*x10d
    for j=2:4*y10d
      dx=lldistkm([lat(p)+latbox(j) lon(p)+lonbox(i-1)],[lat(p)+latbox(j) lon(p)+lonbox(i+1)]).*1000.0; dy=lldistkm([lat(p)+latbox(j-1) lon(p)+lonbox(i)],$
      vor_adv_u(i,j,:) = v_1(i,j,:).*((v_4(i+1,j,:)-v_4(i-1,j,:))./dx);
      vor_adv_v(i,j,:) = v_2(i,j,:).*((v_4(i,j+1,:)-v_4(i,j-1,:))./dy);
    end
  end
  vor_adv=vor_adv_u+vor_adv_v;
  v4=zeros(4*x10d+1,4*y10d+1,length(levs));
  for z=levs
    v4(:,:,z-min(levs)+1)=(vor_adv(:,:,z+1)-vor_adv(:,:,z-1))./((plevs(z+1)-plevs(z-1)).*100.0).*cor(:,:,z-min(levs)+1);
  end
  clear vor_adv_u vor_adv_v vor_adv

% -----------------------------------------------------------------------
% V5: BUOYANCY ADVECTION

  buoyance_adv_u=zeros(4*x10d+1,4*y10d+1,length(plevs)); buoyance_adv_v=zeros(4*x10d+1,4*y10d+1,length(plevs)); buoyance_adv=zeros(4*x10d+1,4*y10d+1,lengt$
  for i=2:4*x10d
    for j=2:4*y10d
      dx=lldistkm([lat(p)+latbox(j) lon(p)+lonbox(i-1)],[lat(p)+latbox(j) lon(p)+lonbox(i+1)]).*1000.0; dy=lldistkm([lat(p)+latbox(j-1) lon(p)+lonbox(i)],$
      buoyance_adv_u(i,j,:) = v_1(i,j,:).*((theta(i+1,j,:) - theta(i-1,j,:))./dx);
      buoyance_adv_v(i,j,:) = v_2(i,j,:).*((theta(i,j+1,:) - theta(i,j-1,:))./dy);
    end
  end
  buoyance_adv=buoyance_adv_u+buoyance_adv_v;
  dbdx=zeros(4*x10d+1,4*y10d+1,length(plevs)); dbdy=zeros(4*x10d+1,4*y10d+1,length(plevs));
  for i = 2:4*x10d
    for j = 2:4*y10d
      dx=lldistkm([lat(p)+latbox(j) lon(p)+lonbox(i-1)],[lat(p)+latbox(j) lon(p)+lonbox(i+1)]).*1000.0; dy=lldistkm([lat(p)+latbox(j-1) lon(p)+lonbox(i)],$
      dbdx(i,j,:) = (buoyance_adv(i+1,j,:)-buoyance_adv(i-1,j,:))./dx;
      dbdy(i,j,:) = (buoyance_adv(i,j+1,:)-buoyance_adv(i,j-1,:))./dy;
    end
  end
  v5=zeros(4*x10d+1,4*y10d+1,length(levs));
  for i = 3:4*x10d-1
    for j = 3:4*y10d-1
      dx=lldistkm([lat(p)+latbox(j) lon(p)+lonbox(i-1)],[lat(p)+latbox(j) lon(p)+lonbox(i+1)]).*1000.0; dy=lldistkm([lat(p)+latbox(j-1) lon(p)+lonbox(i)],$
      v5(i,j,:) = ((dbdx(i+1,j,levs)-dbdx(i-1,j,levs))./dx + (dbdy(i,j+1,levs)-dbdy(i,j-1,levs))./dy).*i_rho(i,j,levs)./theta(i,j,levs);
    end
  end

% -----------------------------------------------------------------------
  v6=v1-v2-v3-v4-v5; % DIABATIC TERM RESIDUAL

% USE A ROTATION MATRIX TO REORIENT EACH VARIABLE RELATIVE TO ITS SHEAR VECTOR.
% THIS IS DESIGNED TO HAVE THE SHEAR VECTOR POINTING UPWARDS ON A PLOT. (downshear --> positive y)
  for z=1:length(levs)
    dwnshr=atan2d(shear_v(p),shear_u(p));
    v1R=imrotate3(v1,dwnshr.*(-1)+90,[0 0 1],'nearest','crop');  % By extracting a larger box initially, we don't accidentally
    v2R=imrotate3(v2,dwnshr.*(-1)+90,[0 0 1],'nearest','crop');  % crop off data that we'll need when considering a 500 km area around the TC.
    v3R=imrotate3(v3,dwnshr.*(-1)+90,[0 0 1],'nearest','crop');
    v4R=imrotate3(v4,dwnshr.*(-1)+90,[0 0 1],'nearest','crop');
    v5R=imrotate3(v5,dwnshr.*(-1)+90,[0 0 1],'nearest','crop');
    v6R=imrotate3(v6,dwnshr.*(-1)+90,[0 0 1],'nearest','crop');
    [x,y]=meshgrid(lon(p)-(2*x10d*xgrid):xgrid:lon(p)+(2*x10d*xgrid),lat(p)-(2*y10d*ygrid):ygrid:lat(p)+(2*y10d*ygrid));
    latR=imrotate(y,dwnshr.*(-1)+90,'nearest','crop');  % Rotate the lats/lons as well, to make sure our distances are correct.
    lonR=imrotate(x,dwnshr.*(-1)+90,'nearest','crop');
    clear x y
    xc=floor(length(xint)/2)+1; yc=floor(length(yint)/2)+1;

% NOW, PLACE VARIABLES INTO BINS BASED ON THE SCHEMES PRESCRIBED INITIALLY
    for t=1:length(tcbins)-1
      for s=1:length(shearbins)-1
        if (shearvar(p) >= shearbins(s) && shearvar(p) < shearbins(s+1) && tcvar(p) >= tcbins(t) && tcvar(p) < tcbins(t+1))
          v1_comp(:,:,z,t,s)=v1_comp(:,:,z,t,s)+v1R(xc-x10d:xc+x10d,yc-y10d:yc+y10d,z);
          v2_comp(:,:,z,t,s)=v2_comp(:,:,z,t,s)+v2R(xc-x10d:xc+x10d,yc-y10d:yc+y10d,z);
          v3_comp(:,:,z,t,s)=v3_comp(:,:,z,t,s)+v3R(xc-x10d:xc+x10d,yc-y10d:yc+y10d,z);
          v4_comp(:,:,z,t,s)=v4_comp(:,:,z,t,s)+v4R(xc-x10d:xc+x10d,yc-y10d:yc+y10d,z);
          v5_comp(:,:,z,t,s)=v5_comp(:,:,z,t,s)+v5R(xc-x10d:xc+x10d,yc-y10d:yc+y10d,z);
          v6_comp(:,:,z,t,s)=v6_comp(:,:,z,t,s)+v6R(xc-x10d:xc+x10d,yc-y10d:yc+y10d,z);
          count(z,t,s)=count(z,t,s)+1; % Accumulate at each vertical level to avoid multiple counting

% NOW THAT WE HAVE THE FULL BOX ACCOUNTED FOR, LET'S CONSIDER THE QUADRANT-SPECIFIC RADIAL PROFILES.
          for q=1:4     % Recall the order: DR --> DL --> UL --> UR
            if (q == 1)
              xrange=xc:xc+x10d; yrange=yc:yc+y10d;   % Extract each corner of the 10-deg shear-rotated box accordingly.
            elseif (q == 2)
              xrange=xc-x10d:xc; yrange=yc:yc+y10d;
            elseif (q == 3)
              xrange=xc-x10d:xc; yrange=yc-y10d:yc;
            elseif (q == 4)
              xrange=xc:xc+x10d; yrange=yc-y10d:yc;
            end
            lat_quad=latR(xrange,yrange); lon_quad=lonR(xrange,yrange);  % Extract lat/lon/data here.
            v1_q=zeros(length(dist),1); v2_q=zeros(length(dist),1); v3_q=zeros(length(dist),1);
            v4_q=zeros(length(dist),1); v5_q=zeros(length(dist),1); v6_q=zeros(length(dist),1);
            distcount=zeros(length(dist),1);
            for x=1:x10d+1
              for y=1:y10d+1     % Assign each point within the quadrant to a radius bin.
                distances=round(lldistkm([lat(p) lon(p)],[lat_quad(x,y) lon_quad(x,y)])/25)*25;  % Rounds to nearest 25 km
                v1_q(find(dist==distances))=v1_q(find(dist==distances))+v1R(xrange(x),yrange(y),z);
                v2_q(find(dist==distances))=v2_q(find(dist==distances))+v2R(xrange(x),yrange(y),z);
                v3_q(find(dist==distances))=v3_q(find(dist==distances))+v3R(xrange(x),yrange(y),z);
                v4_q(find(dist==distances))=v4_q(find(dist==distances))+v4R(xrange(x),yrange(y),z);
                v5_q(find(dist==distances))=v5_q(find(dist==distances))+v5R(xrange(x),yrange(y),z);
                v6_q(find(dist==distances))=v6_q(find(dist==distances))+v6R(xrange(x),yrange(y),z);
                distcount(find(dist==distances))=distcount(find(dist==distances))+1;
              end
            end
            distcount(distcount==0)=NaN;   % If any radii are unaccounted for, change the value to NaN.
            v1_quad(q,:,z,t,s,count(t,s))=v1_q./distcount;  % Take the mean of all points at similar radii within the snapshot
            v2_quad(q,:,z,t,s,count(t,s))=v2_q./distcount;  % but keep the snapshot itself discrete in order to conduct
            v3_quad(q,:,z,t,s,count(t,s))=v3_q./distcount;  % statistical significance testing later on.
            v4_quad(q,:,z,t,s,count(t,s))=v4_q./distcount;
            v5_quad(q,:,z,t,s,count(t,s))=v5_q./distcount;
            v6_quad(q,:,z,t,s,count(t,s))=v6_q./distcount;
          end
          break;
        end
      end
    end
  end
  toc
end

% CALCULATE THE COMPOSITE MEANS OF EACH VARIABLE.
count(count==0)=NaN;
for z=1:length(levs)
  for t=1:length(tcbins)-1
    for s=1:length(shearbins)-1
      v1_comp(:,:,z,t,s)=v1_comp(:,:,z,t,s)./count(z,t,s);
      v2_comp(:,:,z,t,s)=v2_comp(:,:,z,t,s)./count(z,t,s);
      v3_comp(:,:,z,t,s)=v3_comp(:,:,z,t,s)./count(z,t,s);
      v4_comp(:,:,z,t,s)=v4_comp(:,:,z,t,s)./count(z,t,s);
      v5_comp(:,:,z,t,s)=v5_comp(:,:,z,t,s)./count(z,t,s);
      v6_comp(:,:,z,t,s)=v6_comp(:,:,z,t,s)./count(z,t,s);
    end
  end
end

% PLOTTING SECTION - CHANGE PLOT STORAGE DIRECTORY ACCORDINGLY
cd Plots
new_cb=[  0  97 128  % Custom colorbar from nrl_sirkes.rgb in NCL documentation
  0 128 161
  0 161 191
  0 191 224
  0 224 255
  0 255 255
 51 252 252
102 252 252
153 252 252
204 252 252
255 255 255
252 252   0
252 224   0
252 191   0
252 161   0
252 128   0
252  97   0
252  64   0
252  33   0
191   0   0
128   0   0]./255.0;
set(groot,'DefaultFigureColor','white')
set(groot,'DefaultFigureColormap',new_cb)

% PLOT EACH BIN OF SHEAR-RELATIVE VALUES, LOOKING DOWNSHEAR IN THE 10-DEGREE BOX
for t=1:length(tcbins)-1
  for s=1:length(shearbins)-1

%   START WITH SINGLE-LEVEL SPATIAL MAPS
    for z=1:length(levs)
      h=figure;
      imagesc(-5:xgrid:5,-5:ygrid:5,v1_comp(:,:,z,t,s).',[max(max(v1_comp(:,:,z,t,s))).*(-1) max(max(v1_comp(:,:,z,t,s)))]);
      hold on
      for c=1:5  % Plot circles to show radii in 1-deg intervals
        viscircles([0 0],c,'Color','k','LineStyle','--','EnhanceVisibility',false,'LineWidth',1);
      end
      plot([0 0],[-5 5],'--k','LineWidth',1.0);
      plot([-5 5],[0 0],'--k','LineWidth',1.0);
      xlim([-5 5])
      ylim([-5 5])
      set(gca,'XTick',[-5:1:5])
      set(gca,'YTick',[-5:1:5])
      set(gca,'YDir','normal')
      set(gca,'FontSize',14)
      cbh=colorbar;
      caxis([max(max(v1_comp(:,:,z,t,s))).*(-1) max(max(v1_comp(:,:,z,t,s)))]) % May adjust colorbar limits for readability
      xlabel('Across-Shear Degrees from Center')
      ylabel('Along-Shear Degrees from Center')
      str1=strcat(dset," ",num2str(plotlevs(z))," hPa L(\omega)");
      str2=strcat(tcbintext(t),", ",shearbintext(s));
      title([str1,str2])
      xarrow=0; yarrow=0; uarrow=0; varrow=5;
      quiver(xarrow,yarrow,uarrow,varrow,'k','LineWidth',3.5)  % Plot shear vector
      hold off
      pbaspect([1 1 1])
      frameid=sprintf('%02s_LapOmega_Lev_%02d_%02s_%02d_%02s_%02d.png',dset,plotlevs(z),tcplotname,t,shearplotname,s);
      set(gcf,'inverthardcopy','off')
      print(h,'-dpng',frameid);
      clf

      h=figure;
      imagesc(-5:xgrid:5,-5:ygrid:5,v2_comp(:,:,z,t,s).',[max(max(v2_comp(:,:,z,t,s))).*(-1) max(max(v2_comp(:,:,z,t,s)))]);
      hold on
      for c=1:5  % Plot circles to show radii in 1-deg intervals
        viscircles([0 0],c,'Color','k','LineStyle','--','EnhanceVisibility',false,'LineWidth',1);
      end
      plot([0 0],[-5 5],'--k','LineWidth',1.0);
      plot([-5 5],[0 0],'--k','LineWidth',1.0);
      xlim([-5 5])
      ylim([-5 5])
      set(gca,'XTick',[-5:1:5])
      set(gca,'YTick',[-5:1:5])
      set(gca,'YDir','normal')
      set(gca,'FontSize',14)
      cbh=colorbar;
      caxis([max(max(v2_comp(:,:,z,t,s))).*(-1) max(max(v2_comp(:,:,z,t,s)))]) % May adjust colorbar limits for readability
      xlabel('Across-Shear Degrees from Center')
      ylabel('Along-Shear Degrees from Center')
      str1=strcat(dset," ",num2str(plotlevs(z))," hPa DVA_v");
      str2=strcat(tcbintext(t),", ",shearbintext(s));
      title([str1,str2])
      xarrow=0; yarrow=0; uarrow=0; varrow=5;
      quiver(xarrow,yarrow,uarrow,varrow,'k','LineWidth',3.5)  % Plot shear vector
      hold off
      pbaspect([1 1 1])
      frameid=sprintf('%02s_DVAV_Lev_%02d_%02s_%02d_%02s_%02d.png',dset,plotlevs(z),tcplotname,t,shearplotname,s);
      set(gcf,'inverthardcopy','off')
      print(h,'-dpng',frameid);
      clf

      h=figure;
      imagesc(-5:xgrid:5,-5:ygrid:5,v3_comp(:,:,z,t,s).',[max(max(v3_comp(:,:,z,t,s))).*(-1) max(max(v3_comp(:,:,z,t,s)))]);
      hold on
      for c=1:5  % Plot circles to show radii in 1-deg intervals
        viscircles([0 0],c,'Color','k','LineStyle','--','EnhanceVisibility',false,'LineWidth',1);
      end
      plot([0 0],[-5 5],'--k','LineWidth',1.0);
      plot([-5 5],[0 0],'--k','LineWidth',1.0);
      xlim([-5 5])
      ylim([-5 5])
      set(gca,'XTick',[-5:1:5])
      set(gca,'YTick',[-5:1:5])
      set(gca,'YDir','normal')
      set(gca,'FontSize',14)
      cbh=colorbar;
      caxis([max(max(v3_comp(:,:,z,t,s))).*(-1) max(max(v3_comp(:,:,z,t,s)))]) % May adjust colorbar limits for readability
      xlabel('Across-Shear Degrees from Center')
      ylabel('Along-Shear Degrees from Center')
      str1=strcat(dset," ",num2str(plotlevs(z))," hPa Diff. Tilting");
      str2=strcat(tcbintext(t),", ",shearbintext(s));
      title([str1,str2])
      xarrow=0; yarrow=0; uarrow=0; varrow=5;
      quiver(xarrow,yarrow,uarrow,varrow,'k','LineWidth',3.5)  % Plot shear vector
      hold off
      pbaspect([1 1 1])
      frameid=sprintf('%02s_DiffTilt_Lev_%02d_%02s_%02d_%02s_%02d.png',dset,plotlevs(z),tcplotname,t,shearplotname,s);
      set(gcf,'inverthardcopy','off')
      print(h,'-dpng',frameid);
      clf

      h=figure;
      imagesc(-5:xgrid:5,-5:ygrid:5,v4_comp(:,:,z,t,s).',[max(max(v4_comp(:,:,z,t,s))).*(-1) max(max(v4_comp(:,:,z,t,s)))]);
      hold on
      for c=1:5  % Plot circles to show radii in 1-deg intervals
        viscircles([0 0],c,'Color','k','LineStyle','--','EnhanceVisibility',false,'LineWidth',1);
      end
      plot([0 0],[-5 5],'--k','LineWidth',1.0);
      plot([-5 5],[0 0],'--k','LineWidth',1.0);
      xlim([-5 5])
      ylim([-5 5])
      set(gca,'XTick',[-5:1:5])
      set(gca,'YTick',[-5:1:5])
      set(gca,'YDir','normal')
      set(gca,'FontSize',14)
      cbh=colorbar;
      caxis([max(max(v4_comp(:,:,z,t,s))).*(-1) max(max(v4_comp(:,:,z,t,s)))]) % May adjust colorbar limits for readability
      xlabel('Across-Shear Degrees from Center')
      ylabel('Along-Shear Degrees from Center')
      str1=strcat(dset," ",num2str(plotlevs(z))," hPa DVA_h");
      str2=strcat(tcbintext(t),", ",shearbintext(s));
      title([str1,str2])
      xarrow=0; yarrow=0; uarrow=0; varrow=5;
      quiver(xarrow,yarrow,uarrow,varrow,'k','LineWidth',3.5)  % Plot shear vector
      hold off
      pbaspect([1 1 1])
      frameid=sprintf('%02s_DVAH_Lev_%02d_%02s_%02d_%02s_%02d.png',dset,plotlevs(z),tcplotname,t,shearplotname,s);
      set(gcf,'inverthardcopy','off')
      print(h,'-dpng',frameid);
      clf

      h=figure;
      imagesc(-5:xgrid:5,-5:ygrid:5,v5_comp(:,:,z,t,s).',[max(max(v5_comp(:,:,z,t,s))).*(-1) max(max(v5_comp(:,:,z,t,s)))]);
      hold on
      for c=1:5  % Plot circles to show radii in 1-deg intervals
        viscircles([0 0],c,'Color','k','LineStyle','--','EnhanceVisibility',false,'LineWidth',1);
      end
      plot([0 0],[-5 5],'--k','LineWidth',1.0);
      plot([-5 5],[0 0],'--k','LineWidth',1.0);
      xlim([-5 5])
      ylim([-5 5])
      set(gca,'XTick',[-5:1:5])
      set(gca,'YTick',[-5:1:5])
      set(gca,'YDir','normal')
      set(gca,'FontSize',14)
      cbh=colorbar;
      caxis([max(max(v5_comp(:,:,z,t,s))).*(-1) max(max(v5_comp(:,:,z,t,s)))]) % May adjust colorbar limits for readability
      xlabel('Across-Shear Degrees from Center')
      ylabel('Along-Shear Degrees from Center')
      str1=strcat(dset," ",num2str(plotlevs(z))," hPa Buoy. Adv.");
      str2=strcat(tcbintext(t),", ",shearbintext(s));
      title([str1,str2])
      xarrow=0; yarrow=0; uarrow=0; varrow=5;
      quiver(xarrow,yarrow,uarrow,varrow,'k','LineWidth',3.5)  % Plot shear vector
      hold off
      pbaspect([1 1 1])
      frameid=sprintf('%02s_BuoyAdv_Lev_%02d_%02s_%02d_%02s_%02d.png',dset,plotlevs(z),tcplotname,t,shearplotname,s);
      set(gcf,'inverthardcopy','off')
      print(h,'-dpng',frameid);
      clf

      h=figure;
      imagesc(-5:xgrid:5,-5:ygrid:5,v1_comp(:,:,z,t,s).',[max(max(v6_comp(:,:,z,t,s))).*(-1) max(max(v6_comp(:,:,z,t,s)))]);
      hold on
      for c=1:5  % Plot circles to show radii in 1-deg intervals
        viscircles([0 0],c,'Color','k','LineStyle','--','EnhanceVisibility',false,'LineWidth',1);
      end
      plot([0 0],[-5 5],'--k','LineWidth',1.0);
      plot([-5 5],[0 0],'--k','LineWidth',1.0);
      xlim([-5 5])
      ylim([-5 5])
      set(gca,'XTick',[-5:1:5])
      set(gca,'YTick',[-5:1:5])
      set(gca,'YDir','normal')
      set(gca,'FontSize',14)
      cbh=colorbar;
      caxis([max(max(v6_comp(:,:,z,t,s))).*(-1) max(max(v6_comp(:,:,z,t,s)))]) % May adjust colorbar limits for readability
      xlabel('Across-Shear Degrees from Center')
      ylabel('Along-Shear Degrees from Center')
      str1=strcat(dset," ",num2str(plotlevs(z))," hPa Diabatic");
      str2=strcat(tcbintext(t),", ",shearbintext(s));
      title([str1,str2])
      xarrow=0; yarrow=0; uarrow=0; varrow=5;
      quiver(xarrow,yarrow,uarrow,varrow,'k','LineWidth',3.5)  % Plot shear vector
      hold off
      pbaspect([1 1 1])
      frameid=sprintf('%02s_Diab_Lev_%02d_%02s_%02d_%02s_%02d.png',dset,plotlevs(z),tcplotname,t,shearplotname,s);
      set(gcf,'inverthardcopy','off')
      print(h,'-dpng',frameid);
      clf
    end

%   FORMAT AND PLOT QUADRANT-SPECIFIC HEATMAPS HERE. CHOOSE RADIAL RANGE TO AVERAGE OVER.
    for q=1:4
      v1_mean=nanmean(v1_quad(q,:,:,t,s,:),6);
      v2_mean=nanmean(v2_quad(q,:,:,t,s,:),6);
      v3_mean=nanmean(v3_quad(q,:,:,t,s,:),6);
      v4_mean=nanmean(v4_quad(q,:,:,t,s,:),6);
      v5_mean=nanmean(v5_quad(q,:,:,t,s,:),6);
      v6_mean=nanmean(v6_quad(q,:,:,t,s,:),6);
      v_plot=[]; rad=[1:7]; % SET RADIUS RANGE HERE (1-7 --> 0-150 km)
      v_plot(1,:)=nanmean(v1_mean(rad,:),1);
      v_plot(2,:)=nanmean(v2_mean(rad,:),1);
      v_plot(3,:)=nanmean(v3_mean(rad,:),1);
      v_plot(4,:)=nanmean(v4_mean(rad,:),1);
      v_plot(5,:)=nanmean(v5_mean(rad,:),1);
      v_plot(6,:)=nanmean(v6_mean(rad,:),1);

%     KEEPING THE COLOR SCALE LINEAR HERE, BUT THE LOGARITHMIC TRANSFORM APPLIED IN FIG. 9-10 OF
%     PAPER IS AS FOLLOWS: v_trans=sign(v_plot).*(log10(1+abs(v_plot)./(10^(-18))));

      if (q == 1)
        quadrant="Downshear Right Quadrant";
      elseif (q == 2)
        quadrant="Downshear Left Quadrant";
      elseif (q == 3)
        quadrant="Upshear Left Quadrant";
      elseif (q == 4)
        quadrant="Upshear Right Quadrant";
      end
      h=pcolor(0:1:6,plotlevs,v_plot.');
      xlim([0 6])
      ylim([300 950])
      set(gca,'XTick',[0.5:1:5.5])
      set(gca,'XTickLabels',["L(\omega)","DVA_v","DVT","DVA_h","BA","DIAB"])
      set(gca,'YTick',[300:100:900])
      set(gca,'XDir','normal')
      set(gca,'YDir','reverse')
      set(gca,'YScale','log')
      set(gca,'FontSize',14)
      cbh=colorbar;
      caxis([max(max(v_plot)).*(-1) max(max(v_plot))]) % May adjust for readability, just make sure the range is centered on zero for this colorbar!
      grid on
      xlabel('Process')
      ylabel('Pressure (hPa)')
      str1=strcat(dset," \omega Equation - ",num2str(dist(rad(1))),"-",num2str(dist(rad(end)))," km");
      str2=strcat(tcbintext(t),", ",shearbintext(s));
      str3=quadrant;
      title([str1,str2,str3])
      pbaspect([2 1 1])   % Makes the plot taller than it is wide, so multiple quadrant-specific profiles can fit on one page legibly.
      frameid=sprintf('%02s_OmegaEquation_InnerCore_%02s_%02d_%02s_%02d_%02d.png',dset,tcplotname,t,shearplotname,s,q);  % Change plot title accordingly with radius range setting.
      set(gcf,'inverthardcopy','off')
      print(gcf,'-dpng',frameid);
      clf

%     PLOT RESULTS AS LINES
      h=figure;
      plot(v_plot(1,:),plotlevs,'Color','k','LineStyle','-','LineWidth',1.8)
      hold on
      plot(v_plot(2,:),plotlevs,'Color',[0,1,0],'LineStyle','-','LineWidth',1.8)
      plot(v_plot(3,:),plotlevs,'Color',[0,1,1],'LineStyle','-','LineWidth',1.8)
      plot(v_plot(4,:),plotlevs,'Color',[1,0.5,0],'LineStyle','-','LineWidth',1.8)
      plot(v_plot(5,:),plotlevs,'Color',[0,0,1],'LineStyle','-','LineWidth',1.8)
      plot(v_plot(6,:),plotlevs,'Color',[1,0,1],'LineStyle','-','LineWidth',1.8)
      ylim([300 950])
      ylabel('Pressure (hPa)')
      set(gca,'YTick',[300:100:900])
      xlabel('Contribution to \omega (Pa^{-1} s^{-3})')
      leg=legend('L(\omega)','DVA_v','DVT','DVA_h','BA','DIAB','location','east');
      leg.FontSize=13;
      set(gca,'XDir','normal')
      set(gca,'YDir','reverse')
      set(gca,'YScale','log')
      set(gca,'FontSize',16)
      grid on
      str1=strcat(dset," \omega Equation - ",num2str(dist(rad(1))),"-",num2str(dist(rad(end)))," km");
      str2=strcat(tcbintext(t),", ",shearbintext(s));
      str3=quadrant;
      title([str1,str2,str3])
      pbaspect([2 1 1])
      frameid=sprintf('%02s_OmegaEquation_InnerCore_LINES_%02d.png',dset,q);
      set(gcf,'inverthardcopy','off')
      print(gcf,'-dpng',frameid);
      clf
    end
  end
end
