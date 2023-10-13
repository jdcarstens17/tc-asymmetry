% FUNCTION Thermo - Edited by Jake Carstens 10/12/2023 for manuscript on tropical cyclone asymmetry in reanalyses
%                   Generalized version of script used to generate Figure 6 of manuscript.
% PURPOSE: Compute and plot spatial maps of specific/relative humidity, vertical motion, surface heat fluxes,
%          and thermal anomalies in shear-relative space, at the vertical level of the user's choice.
% NOTES: 1. Storm-centered data snapshots, TempestExtremes-derived TC tracks, and wind shear files will be provided
%           in same directory within repository, so filepaths should not need to be altered.
%        2. The user has the choice to use either .mat or .nc format for the wind shear and storm-centered data. TempestExtremes
%           tracks are provided as .txt files.
% PROCEDURE:
% 1. Choose reanalysis dataset (currently ERA5 or CFSR), loading in its TC tracks, associated wind shear,
%    and information about the grid.
% 2. Pre-allocate matrices for specific humidity (v1), equivalent potential temperature (v2), potential
%    anomaly (v3), vertical velocity (v4), layer-average relative humidity (v5), latent heat flux (v6),
%    sensible heat flux (v7), binned by both TC behavior and wind shear magnitude.
% 3. Loop through all TC snapshots, rotate data to shear vector, and place into TC/shear bin composite at each vertical level.
% 4. Plot each variable as composite spatial maps. In general, humidity/surface fluxes (and its azimuthal anomalies)
%    will be shaded, while temperature azimuthal anomalies and vertical velocity will be contoured.

% FIRST, CHOOSE THE DATASET. THIS INCLUDES INFORMATION ABOUT THE GRID SPACING, USEFUL FOR DRAWING A 10-DEGREE BOX
% AROUND THE TC (x10d/y10d), AS WELL AS RE-CENTERING (xtotal) THE DOMAIN ABOUT THE TC.
dset="CFSR";
if (dset == "CFSR")
  x10d=10; y10d=10; xtotal=720; xgrid=0.50; ygrid=0.50; dset_c='CFSR';
elseif (dset == "ERA5")
  x10d=20; y10d=20; xtotal=1440; xgrid=0.25; ygrid=0.25; dset_c='ERA5';
end
lev=[100:25:250 300:50:750 775:25:1000];   % Vertical levels, which are the same in both reanalyses here. 27 total.

% CHOOSE VERTICAL LEVELS TO PLOT DATA FOR, INCLUDING LAYER BOUNDARIES FOR RELATIVE HUMIDITY CALCULATION.
% CURRENT SETTINGS: Q - 925 hPa; theta-E - 925 hPa; theta - 600 hPa; w - 600 hPa; RH - 700-500 hPa
lev_sh=24; lev_thetae=24; lev_theta=14; lev_w=14; lev_rh=12:16;

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
% v1: Specific humidity  v2: Theta-E anomaly  v3: Theta anomaly  v4: Omega
% v5: Relative humidity  v6: Latent heat flux v7: Sensible heat flux
count=zeros(length(tcbins)-1,length(shearbins)-1);                      % Accumulate all snapshots in a particular TC/shear bin.
v1_comp=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);  % Draw a 10-deg by 10-deg box for the TC composite.
v2_comp=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);
v3_comp=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);
v4_comp=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);
v5_comp=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);
v6_comp=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);
v7_comp=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);

% IF CHOOSING THE .mat OPTION, QUEUE UP THE FILES HERE
W=matfile(['StormCenteredData/' dset_c '_W.mat']);
T=matfile(['StormCenteredData/' dset_c '_T.mat']);
Q=matfile(['StormCenteredData/' dset_c '_Q.mat']);
LHF=matfile(['StormCenteredData/' dset_c '_LHF.mat']);
SHF=matfile(['StormCenteredData/' dset_c '_SHF.mat']);

% ENTER LOOP MOVING THROUGH ALL TC SNAPSHOTS, WHERE RELEVANT DATA WILL BE LOADED, ROTATED ABOUT THE SHEAR, AND PLACED INTO COMPOSITES.
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
  sh=Q.sh(:,:,:,p);     % Specific humidity
  temp=T.temp(:,:,:,p); % Temperature
  w=W.w(:,:,:,p);       % Vertical velocity
  lhf=LHF.lhf(:,:,p);   % Latent heat flux
  shf=SHF.shf(:,:,p);   % Sensible heat flux

% .nc SECTION
  start=[1 1 1 p]; interval=[length([-10:xgrid:10]) length([-10:ygrid:10]) length(lev) 1];
  sh=ncread(['StormCenteredData/' dset_c '_Q.nc'],'sh',start,interval);
  temp=ncread(['StormCenteredData/' dset_c '_T.nc'],'temp',start,interval);
  w=ncread(['StormCenteredData/' dset_c '_W.nc'],'w',start,interval);
  lhf=ncread(['StormCenteredData/' dset_c '_LHF.nc'],'lhf',start,interval);
  shf=ncread(['StormCenteredData/' dset_c '_SHF.nc'],'shf',start,interval);

  v1=sh(:,:,lev_sh).*1000.0;                                % SPECIFIC HUMIDITY
  mix=sh(:,:,lev_sh)./(1-sh(:,:,lev_sh)); L=2.5e6; cp=1005; % THETA-E
  v2=(temp(:,:,lev_thetae)+(L./cp.*mix)).*((1000./lev(lev_thetae)).^0.286);
  v3=temp(:,:,lev_theta).*((1000./lev(lev_theta)).^0.286);  % THETA
  v4=w(:,:,lev_w);                                          % OMEGA

  v5=zeros(4*x10d+1,4*y10d+1);                              % RELATIVE HUMIDITY
  for z=1:length(lev_rh)
    nu=0.263.*lev(lev_rh(z)).*sh(:,:,lev_rh(z));
    de=exp(17.67.*(temp(:,:,lev_rh(z))-273.15)./(273.15-29.65));
    v5=v5+(nu./de);
  end
  v5=v5./length(lev_rh).*100.0;
  if (dset == "ERA5")                                       % LATENT AND SENSIBLE HEAT FLUX
    v6=lhf(-1); v7=shf.*(-1);
  elseif (dset == "CFSR")
    v6=lhf; v7=shf;
  end

% USE A ROTATION MATRIX TO REORIENT EACH VARIABLE RELATIVE TO ITS SHEAR VECTOR.
% THIS IS DESIGNED TO HAVE THE SHEAR VECTOR POINTING UPWARDS ON A PLOT. (downshear --> positive y)
  dwnshr=atan2d(shear_v(p),shear_u(p));
  v1R=imrotate(v1,dwnshr.*(-1)+90,'nearest','crop');  % By extracting a larger box initially, we don't accidentally
  v2R=imrotate(v2,dwnshr.*(-1)+90,'nearest','crop');  % crop off data that we'll need when considering a 500 km area around the TC.
  v3R=imrotate(v3,dwnshr.*(-1)+90,'nearest','crop');
  v4R=imrotate(v4,dwnshr.*(-1)+90,'nearest','crop');
  v5R=imrotate(v5,dwnshr.*(-1)+90,'nearest','crop');
  v6R=imrotate(v6,dwnshr.*(-1)+90,'nearest','crop');
  v7R=imrotate(v7,dwnshr.*(-1)+90,'nearest','crop');
  [x,y]=meshgrid(lon(p)-(2*x10d*xgrid):xgrid:lon(p)+(2*x10d*xgrid),lat(p)-(2*y10d*ygrid):ygrid:lat(p)+(2*y10d*ygrid));
  latR=imrotate(y,dwnshr.*(-1)+90,'nearest','crop');  % Rotate the lats/lons as well, to make sure our distances are correct.
  lonR=imrotate(x,dwnshr.*(-1)+90,'nearest','crop');
  clear x y
  xc=floor((4*x10d+1)/2)+1; yc=floor((4*y10d+1)/2)+1; % Extract the center of the new rotated square as the TC center to separate quadrants.

% NOW, PLACE VARIABLES INTO BINS BASED ON THE SCHEMES PRESCRIBED INITIALLY
  for t=1:length(tcbins)-1
    for s=1:length(shearbins)-1
      if (shearvar(p) >= shearbins(s) && shearvar(p) < shearbins(s+1) && tcvar(p) >= tcbins(t) && tcvar(p) < tcbins(t+1))
        v1_comp(:,:,t,s)=v1_comp(:,:,t,s)+v1R(xc-x10d:xc+x10d,yc-y10d:yc+y10d);  % Now just extract the 10-deg shear-rotated box.
        v2_comp(:,:,t,s)=v2_comp(:,:,t,s)+v2R(xc-x10d:xc+x10d,yc-y10d:yc+y10d);
        v3_comp(:,:,t,s)=v3_comp(:,:,t,s)+v3R(xc-x10d:xc+x10d,yc-y10d:yc+y10d);
        v4_comp(:,:,t,s)=v4_comp(:,:,t,s)+v4R(xc-x10d:xc+x10d,yc-y10d:yc+y10d);
        v5_comp(:,:,t,s)=v5_comp(:,:,t,s)+v5R(xc-x10d:xc+x10d,yc-y10d:yc+y10d);
        v6_comp(:,:,t,s)=v6_comp(:,:,t,s)+v6R(xc-x10d:xc+x10d,yc-y10d:yc+y10d);
        v7_comp(:,:,t,s)=v7_comp(:,:,t,s)+v7R(xc-x10d:xc+x10d,yc-y10d:yc+y10d);
        count(t,s)=count(t,s)+1;        % Accumulate each snapshot falling into the bin, since we'll calculate a composite mean later.
      end
    end
  end
  toc
end

% CALCULATE THE COMPOSITE MEANS OF EACH VARIABLE.
count(count==0)=NaN;
for t=1:length(tcbins)-1
  for s=1:length(shearbins)-1
    v1_comp(:,:,t,s)=v1_comp(:,:,t,s)./count(t,s);
    v2_comp(:,:,t,s)=v2_comp(:,:,t,s)./count(t,s);
    v3_comp(:,:,t,s)=v3_comp(:,:,t,s)./count(t,s);
    v4_comp(:,:,t,s)=v4_comp(:,:,t,s)./count(t,s);
    v5_comp(:,:,t,s)=v5_comp(:,:,t,s)./count(t,s);
    v6_comp(:,:,t,s)=v6_comp(:,:,t,s)./count(t,s);
    v7_comp(:,:,t,s)=v7_comp(:,:,t,s)./count(t,s);
  end
end

% SECTION TO REWRITE SURFACE FLUXES, SPECIFIC/RELATIVE HUMIDITY, AND THETA/THETA-E AS ANOMALIES FROM THE AZIMUTHAL MEAN

v1_new=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);
v2_new=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);
v3_new=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);
v5_new=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);
v6_new=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);
v7_new=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);
distancesnew=zeros(2*x10d+1,2*y10d+1); xc=floor((2*x10d+1)/2)+1; yc=floor((2*y10d+1)/2)+1;
for a=1:2*x10d+1
  for b=1:2*y10d+1
    distancesnew(a,b)=sqrt((xc-a).^2+(yc-b).^2);
  end
end
for t=1:length(tcbins)-1
  for s=1:length(shearbins)-1
    v1_bin=v1_comp(:,:,t,s);
    v2_bin=v2_comp(:,:,t,s);
    v3_bin=v3_comp(:,:,t,s);
    v5_bin=v5_comp(:,:,t,s);
    v6_bin=v6_comp(:,:,t,s);
    v7_bin=v7_comp(:,:,t,s);
    for a=1:2*x10d+1
      for b=1:2*y10d+1
        distnew=sqrt((xc-a).^2+(yc-b).^2);
        v1_new(a,b,t,s)=v1_bin(a,b)-mean(mean(v1_bin(find(distancesnew==distnew))));
        v2_new(a,b,t,s)=v2_bin(a,b)-mean(mean(v2_bin(find(distancesnew==distnew))));
        v3_new(a,b,t,s)=v3_bin(a,b)-mean(mean(v3_bin(find(distancesnew==distnew))));
        v5_new(a,b,t,s)=v5_bin(a,b)-mean(mean(v5_bin(find(distancesnew==distnew))));
        v6_new(a,b,t,s)=v6_bin(a,b)-mean(mean(v6_bin(find(distancesnew==distnew))));
        v7_new(a,b,t,s)=v7_bin(a,b)-mean(mean(v7_bin(find(distancesnew==distnew))));
      end
    end
  end
end

% PLOTTING SECTION - CHANGE PLOT STORAGE DIRECTORY ACCORDINGLY
cd Plots
set(groot,'DefaultFigureColor','white')

% EDIT BELOW THIS POINT THURSDAY EVENING!!!!!!!!!

% PLOT EACH BIN OF SHEAR-RELATIVE VALUES, LOOKING DOWNSHEAR IN THE 10-DEGREE BOX
for t=1:length(tcbins)-1
  for s=1:length(shearbins)-1

  % START BY PLOTTING SPECIFIC HUMIDITY WITH THETA-E ANOMALY OVERLAID
    set(groot,'DefaultFigureColormap',flipud(parula))
    h=figure;
    imagesc(-5:xgrid:5,-5:ygrid:5,v1_comp(:,:,t,s).',[min(min(v1_comp(:,:,t,s))) max(max(v1_comp(:,:,t,s)))]);
    hold on
    for c=1:5   % Plot circles to better visualize radial distances from center.
      viscircles([0 0],c,'Color','k','LineStyle','--','EnhanceVisibility',false,'LineWidth',1);
    end
    plot([0 0],[-5 5],'--k','LineWidth',1.0);
    plot([-5 5],[0 0],'--k','LineWidth',1.0);  % Theta-e anomalies plotted in 0.5 K intervals, may change depending on preference/level.
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[-2.5 -2.5],'LineColor','b','LineStyle','-','LineWidth',2.5)
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[-2 -2],'LineColor','b','LineStyle','-','LineWidth',2.2)
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[-1.5 -1.5],'LineColor','b','LineStyle','-','LineWidth',1.9)
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[-1 -1],'LineColor','b','LineStyle','-','LineWidth',1.6)
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[-0.5 -0.5],'LineColor','b','LineStyle','-','LineWidth',1.3)
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[0.5 0.5],'LineColor','r','LineStyle','-','LineWidth',1.3)
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[1 1],'LineColor','r','LineStyle','-','LineWidth',1.6)
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[1.5 1.5],'LineColor','r','LineStyle','-','LineWidth',1.9)
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[2 2],'LineColor','r','LineStyle','-','LineWidth',2.2)
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[2.5 2.5],'LineColor','r','LineStyle','-','LineWidth',2.5)  % Thicker --> greater anomaly
    xlim([-5 5])
    ylim([-5 5])
    set(gca,'XTick',[-5:1:5])
    set(gca,'YTick',[-5:1:5])
    set(gca,'YDir','normal')
    set(gca,'FontSize',14)
    cbh=colorbar;
    caxis([min(min(v1_comp(:,:,t,s))) max(max(v1_comp(:,:,t,s)))])   % May adjust colorbar limits for readability
    xlabel('Across-Shear Degrees from Center')
    ylabel('Along-Shear Degrees from Center')
    str1=strcat(dset," ",num2str(lev(lev_sh))," hPa Specific Humidity (g/kg)");  % All necessary info will be included in plot title, and name below
    str2=strcat(tcbintext(t),", ",shearbintext(s));
    title([str1,str2])
    xarrow=0; yarrow=0; uarrow=0; varrow=5;
    quiver(xarrow,yarrow,uarrow,varrow,'k','LineWidth',3.5)  % Plot shear vector
    hold off
    pbaspect([1 1 1])
    frameid=sprintf('%02s_Humidity_Lev%02d_%02s_%02d_%02s_%02d.png',dset,lev_sh,tcplotname,t,shearplotname,s);
    set(gcf,'inverthardcopy','off')
    print(h,'-dpng',frameid);
    clf

  % SAME PLOT AS ABOVE, BUT PLOTTING THE SPECIFIC HUMIDITY AZIMUTHAL ANOMALY INSTEAD
    set(groot,'DefaultFigureColormap',flipud(parula))
    h=figure;
    imagesc(-5:xgrid:5,-5:ygrid:5,v1_new(:,:,t,s).',[min(min(v1_new(:,:,t,s))) max(max(v1_new(:,:,t,s)))]);
    hold on
    for c=1:5   % Plot circles to better visualize radial distances from center.
      viscircles([0 0],c,'Color','k','LineStyle','--','EnhanceVisibility',false,'LineWidth',1);
    end
    plot([0 0],[-5 5],'--k','LineWidth',1.0);
    plot([-5 5],[0 0],'--k','LineWidth',1.0);  % Theta-e anomalies plotted in 0.5 K intervals, may change depending on preference/level.
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[-2.5 -2.5],'LineColor','b','LineStyle','-','LineWidth',2.5)
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[-2 -2],'LineColor','b','LineStyle','-','LineWidth',2.2)
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[-1.5 -1.5],'LineColor','b','LineStyle','-','LineWidth',1.9)
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[-1 -1],'LineColor','b','LineStyle','-','LineWidth',1.6)
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[-0.5 -0.5],'LineColor','b','LineStyle','-','LineWidth',1.3)
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[0.5 0.5],'LineColor','r','LineStyle','-','LineWidth',1.3)
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[1 1],'LineColor','r','LineStyle','-','LineWidth',1.6)
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[1.5 1.5],'LineColor','r','LineStyle','-','LineWidth',1.9)
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[2 2],'LineColor','r','LineStyle','-','LineWidth',2.2)
    contour(-5:xgrid:5,-5:ygrid:5,v2_new(:,:,t,s).',[2.5 2.5],'LineColor','r','LineStyle','-','LineWidth',2.5)  % Thicker --> greater anomaly
    xlim([-5 5])
    ylim([-5 5])
    set(gca,'XTick',[-5:1:5])
    set(gca,'YTick',[-5:1:5])
    set(gca,'YDir','normal')
    set(gca,'FontSize',14)
    cbh=colorbar;
    caxis([min(min(v1_comp(:,:,t,s))) max(max(v1_comp(:,:,t,s)))])   % May adjust colorbar limits for readability
    xlabel('Across-Shear Degrees from Center')
    ylabel('Along-Shear Degrees from Center')
    str1=strcat(dset," ",num2str(lev(lev_sh))," hPa Specific Humidity Anomaly (g/kg)");  % All necessary info will be included in plot title, and name below
    str2=strcat(tcbintext(t),", ",shearbintext(s));
    title([str1,str2])
    xarrow=0; yarrow=0; uarrow=0; varrow=5;
    quiver(xarrow,yarrow,uarrow,varrow,'k','LineWidth',3.5)  % Plot shear vector
    hold off
    pbaspect([1 1 1])
    frameid=sprintf('%02s_HumidityAnom_Lev%02d_%02s_%02d_%02s_%02d.png',dset,lev_sh,tcplotname,t,shearplotname,s);
    set(gcf,'inverthardcopy','off')
    print(h,'-dpng',frameid);
    clf

% NEXT: VERTICAL VELOCITY WITH THETA ANOMALY OVERLAID
    set(groot,'DefaultFigureColormap',flipud(hot))
    h=figure;
    imagesc(-5:xgrid:5,-5:ygrid:5,v4_comp(:,:,t,s).',[min(min(v4_comp(:,:,t,s))) 0]);
    hold on
    for c=1:5   % Plot circles to better visualize radial distances from center.
      viscircles([0 0],c,'Color','k','LineStyle','--','EnhanceVisibility',false,'LineWidth',1);
    end
    plot([0 0],[-5 5],'--k','LineWidth',1.0);
    plot([-5 5],[0 0],'--k','LineWidth',1.0);  % Theta anomalies plotted in 0.2 K intervals
    contour(-5:xgrid:5,-5:ygrid:5,v3_new(:,:,t,s).',[-1 -1],'LineColor','b','LineStyle','-','LineWidth',2.5)
    contour(-5:xgrid:5,-5:ygrid:5,v3_new(:,:,t,s).',[-0.8 -0.8],'LineColor','b','LineStyle','-','LineWidth',2.2)
    contour(-5:xgrid:5,-5:ygrid:5,v3_new(:,:,t,s).',[-0.6 -0.6],'LineColor','b','LineStyle','-','LineWidth',1.9)
    contour(-5:xgrid:5,-5:ygrid:5,v3_new(:,:,t,s).',[-0.4 -0.4],'LineColor','b','LineStyle','-','LineWidth',1.6)
    contour(-5:xgrid:5,-5:ygrid:5,v3_new(:,:,t,s).',[-0.2 -0.2],'LineColor','b','LineStyle','-','LineWidth',1.3)
    contour(-5:xgrid:5,-5:ygrid:5,v3_new(:,:,t,s).',[0.2 0.2],'LineColor','r','LineStyle','-','LineWidth',1.3)
    contour(-5:xgrid:5,-5:ygrid:5,v3_new(:,:,t,s).',[0.4 0.4],'LineColor','r','LineStyle','-','LineWidth',1.6)
    contour(-5:xgrid:5,-5:ygrid:5,v3_new(:,:,t,s).',[0.6 0.6],'LineColor','r','LineStyle','-','LineWidth',1.9)
    contour(-5:xgrid:5,-5:ygrid:5,v3_new(:,:,t,s).',[0.8 0.8],'LineColor','r','LineStyle','-','LineWidth',2.2)
    contour(-5:xgrid:5,-5:ygrid:5,v3_new(:,:,t,s).',[1 1],'LineColor','r','LineStyle','-','LineWidth',2.5)  % Thicker --> greater anomaly
    xlim([-5 5])
    ylim([-5 5])
    set(gca,'XTick',[-5:1:5])
    set(gca,'YTick',[-5:1:5])
    set(gca,'YDir','normal')
    set(gca,'FontSize',14)
    cbh=colorbar;
    caxis([min(min(v4_comp(:,:,t,s))) 0])   % May adjust colorbar limits for readability
    xlabel('Across-Shear Degrees from Center')
    ylabel('Along-Shear Degrees from Center')
    str1=strcat(dset," ",num2str(lev(lev_w))," hPa \omega (Pa/s)");  % All necessary info will be included in plot title, and name below
    str2=strcat(tcbintext(t),", ",shearbintext(s));
    title([str1,str2])
    xarrow=0; yarrow=0; uarrow=0; varrow=5;
    quiver(xarrow,yarrow,uarrow,varrow,'k','LineWidth',3.5)  % Plot shear vector
    hold off
    pbaspect([1 1 1])
    frameid=sprintf('%02s_Omega_Lev%02d_%02s_%02d_%02s_%02d.png',dset,lev_w,tcplotname,t,shearplotname,s);
    set(gcf,'inverthardcopy','off')
    print(h,'-dpng',frameid);
    clf

% NEXT: LATENT HEAT FLUX WITH LHF ANOMALY OVERLAID
    set(groot,'DefaultFigureColormap',flipud(hot))
    h=figure;
    imagesc(-5:xgrid:5,-5:ygrid:5,v6_comp(:,:,t,s).',[min(min(v6_comp(:,:,t,s))) max(max(v6_comp(:,:,t,s)))]);
    hold on
    for c=1:5   % Plot circles to better visualize radial distances from center.
      viscircles([0 0],c,'Color','k','LineStyle','--','EnhanceVisibility',false,'LineWidth',1);
    end
    plot([0 0],[-5 5],'--k','LineWidth',1.0);
    plot([-5 5],[0 0],'--k','LineWidth',1.0);  % LHF anomalies plotted in 10 W/m^2 intervals
    contour(-5:xgrid:5,-5:ygrid:5,v6_new(:,:,t,s).',[-50 -50],'LineColor','b','LineStyle','-','LineWidth',2.5)
    contour(-5:xgrid:5,-5:ygrid:5,v6_new(:,:,t,s).',[-40 -40],'LineColor','b','LineStyle','-','LineWidth',2.2)
    contour(-5:xgrid:5,-5:ygrid:5,v6_new(:,:,t,s).',[-30 -30],'LineColor','b','LineStyle','-','LineWidth',1.9)
    contour(-5:xgrid:5,-5:ygrid:5,v6_new(:,:,t,s).',[-20 -20],'LineColor','b','LineStyle','-','LineWidth',1.6)
    contour(-5:xgrid:5,-5:ygrid:5,v6_new(:,:,t,s).',[-10 -10],'LineColor','b','LineStyle','-','LineWidth',1.3)
    contour(-5:xgrid:5,-5:ygrid:5,v6_new(:,:,t,s).',[10 10],'LineColor','r','LineStyle','-','LineWidth',1.3)
    contour(-5:xgrid:5,-5:ygrid:5,v6_new(:,:,t,s).',[20 20],'LineColor','r','LineStyle','-','LineWidth',1.6)
    contour(-5:xgrid:5,-5:ygrid:5,v6_new(:,:,t,s).',[30 30],'LineColor','r','LineStyle','-','LineWidth',1.9)
    contour(-5:xgrid:5,-5:ygrid:5,v6_new(:,:,t,s).',[40 40],'LineColor','r','LineStyle','-','LineWidth',2.2)
    contour(-5:xgrid:5,-5:ygrid:5,v6_new(:,:,t,s).',[50 50],'LineColor','r','LineStyle','-','LineWidth',2.5)  % Thicker --> greater anomaly
    xlim([-5 5])
    ylim([-5 5])
    set(gca,'XTick',[-5:1:5])
    set(gca,'YTick',[-5:1:5])
    set(gca,'YDir','normal')
    set(gca,'FontSize',14)
    cbh=colorbar;
    caxis([min(min(v6_comp(:,:,t,s))) max(max(v6_comp(:,:,t,s)))])   % May adjust colorbar limits for readability
    xlabel('Across-Shear Degrees from Center')
    ylabel('Along-Shear Degrees from Center')
    str1=strcat(dset," Latent Heat Flux (W/m^2)");  % All necessary info will be included in plot title, and name below
    str2=strcat(tcbintext(t),", ",shearbintext(s));
    title([str1,str2])
    xarrow=0; yarrow=0; uarrow=0; varrow=5;
    quiver(xarrow,yarrow,uarrow,varrow,'k','LineWidth',3.5)  % Plot shear vector
    hold off
    pbaspect([1 1 1])
    frameid=sprintf('%02s_LHF_%02s_%02d_%02s_%02d.png',dset,tcplotname,t,shearplotname,s);
    set(gcf,'inverthardcopy','off')
    print(h,'-dpng',frameid);
    clf

% NEXT: LATENT HEAT FLUX WITH LHF ANOMALY OVERLAID
    set(groot,'DefaultFigureColormap',flipud(hot))
    h=figure;
    imagesc(-5:xgrid:5,-5:ygrid:5,v7_comp(:,:,t,s).',[min(min(v7_comp(:,:,t,s))) max(max(v7_comp(:,:,t,s)))]);
    hold on
    for c=1:5   % Plot circles to better visualize radial distances from center.
      viscircles([0 0],c,'Color','k','LineStyle','--','EnhanceVisibility',false,'LineWidth',1);
    end
    plot([0 0],[-5 5],'--k','LineWidth',1.0);
    plot([-5 5],[0 0],'--k','LineWidth',1.0);  % SHF anomalies plotted in 10 W/m^2 intervals
    contour(-5:xgrid:5,-5:ygrid:5,v7_new(:,:,t,s).',[-50 -50],'LineColor','b','LineStyle','-','LineWidth',2.5)
    contour(-5:xgrid:5,-5:ygrid:5,v7_new(:,:,t,s).',[-40 -40],'LineColor','b','LineStyle','-','LineWidth',2.2)
    contour(-5:xgrid:5,-5:ygrid:5,v7_new(:,:,t,s).',[-30 -30],'LineColor','b','LineStyle','-','LineWidth',1.9)
    contour(-5:xgrid:5,-5:ygrid:5,v7_new(:,:,t,s).',[-20 -20],'LineColor','b','LineStyle','-','LineWidth',1.6)
    contour(-5:xgrid:5,-5:ygrid:5,v7_new(:,:,t,s).',[-10 -10],'LineColor','b','LineStyle','-','LineWidth',1.3)
    contour(-5:xgrid:5,-5:ygrid:5,v7_new(:,:,t,s).',[10 10],'LineColor','r','LineStyle','-','LineWidth',1.3)
    contour(-5:xgrid:5,-5:ygrid:5,v7_new(:,:,t,s).',[20 20],'LineColor','r','LineStyle','-','LineWidth',1.6)
    contour(-5:xgrid:5,-5:ygrid:5,v7_new(:,:,t,s).',[30 30],'LineColor','r','LineStyle','-','LineWidth',1.9)
    contour(-5:xgrid:5,-5:ygrid:5,v7_new(:,:,t,s).',[40 40],'LineColor','r','LineStyle','-','LineWidth',2.2)
    contour(-5:xgrid:5,-5:ygrid:5,v7_new(:,:,t,s).',[50 50],'LineColor','r','LineStyle','-','LineWidth',2.5)  % Thicker --> greater anomaly
    xlim([-5 5])
    ylim([-5 5])
    set(gca,'XTick',[-5:1:5])
    set(gca,'YTick',[-5:1:5])
    set(gca,'YDir','normal')
    set(gca,'FontSize',14)
    cbh=colorbar;
    caxis([min(min(v7_comp(:,:,t,s))) max(max(v7_comp(:,:,t,s)))])   % May adjust colorbar limits for readability
    xlabel('Across-Shear Degrees from Center')
    ylabel('Along-Shear Degrees from Center')
    str1=strcat(dset,"Sensible Heat Flux (W/m^2)");  % All necessary info will be included in plot title, and name below
    str2=strcat(tcbintext(t),", ",shearbintext(s));
    title([str1,str2])
    xarrow=0; yarrow=0; uarrow=0; varrow=5;
    quiver(xarrow,yarrow,uarrow,varrow,'k','LineWidth',3.5)  % Plot shear vector
    hold off
    pbaspect([1 1 1])
    frameid=sprintf('%02s_SHF_%02s_%02d_%02s_%02d.png',dset,tcplotname,t,shearplotname,s);
    set(gcf,'inverthardcopy','off')
    print(h,'-dpng',frameid);
    clf

% NEXT: RELATIVE HUMIDITY WITH THETA AND VERTICAL VELOCITY OVERLAID
    set(groot,'DefaultFigureColormap',flipud(parula))
    h=figure;
    imagesc(-5:xgrid:5,-5:ygrid:5,v5_comp(:,:,t,s).',[min(min(v5_comp(:,:,t,s))) max(max(v5_comp(:,:,t,s)))]);
    hold on
    for c=1:5   % Plot circles to better visualize radial distances from center.
      viscircles([0 0],c,'Color','k','LineStyle','--','EnhanceVisibility',false,'LineWidth',1);
    end
    plot([0 0],[-5 5],'--k','LineWidth',1.0);
    plot([-5 5],[0 0],'--k','LineWidth',1.0);  % W in 0.5 Pa/s intervals (ascent only); theta anomalies in 0.2 K intervals
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-3 -3],'LineColor','k','LineStyle','-','LineWidth',2.5)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-2.5 -2.5],'LineColor','k','LineStyle','-','LineWidth',2.2)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-2 -2],'LineColor','k','LineStyle','-','LineWidth',1.9)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-1.5 -1.5],'LineColor','k','LineStyle','-','LineWidth',1.6)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-1 -1],'LineColor','k','LineStyle','-','LineWidth',1.3)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-0.5 -0.5],'LineColor','k','LineStyle','-','LineWidth',1.0)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[-1 -1],'LineColor','b','LineStyle','-','LineWidth',2.5)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[-0.8 -0.8],'LineColor','b','LineStyle','-','LineWidth',2.2)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[-0.6 -0.6],'LineColor','b','LineStyle','-','LineWidth',1.9)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[-0.4 -0.4],'LineColor','b','LineStyle','-','LineWidth',1.6)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[-0.2 -0.2],'LineColor','b','LineStyle','-','LineWidth',1.3)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[0.2 0.2],'LineColor','r','LineStyle','-','LineWidth',1.3)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[0.4 0.4],'LineColor','r','LineStyle','-','LineWidth',1.6)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[0.6 0.6],'LineColor','r','LineStyle','-','LineWidth',1.9)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[0.8 0.8],'LineColor','r','LineStyle','-','LineWidth',2.2)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[1 1],'LineColor','r','LineStyle','-','LineWidth',2.5)
    xlim([-5 5])
    ylim([-5 5])
    set(gca,'XTick',[-5:1:5])
    set(gca,'YTick',[-5:1:5])
    set(gca,'YDir','normal')
    set(gca,'FontSize',14)
    cbh=colorbar;
    caxis([min(min(v5_comp(:,:,t,s))) max(max(v5_comp(:,:,t,s)))])   % May adjust colorbar limits for readability
    xlabel('Across-Shear Degrees from Center')
    ylabel('Along-Shear Degrees from Center')
    str1=strcat(dset," ",lev(lev_rh(end)),"-",lev(lev_rh(1))," hPa Relative Humidity (%)");  % All necessary info will be included in plot title, and name below
    str2=strcat(tcbintext(t),", ",shearbintext(s));
    title([str1,str2])
    xarrow=0; yarrow=0; uarrow=0; varrow=5;
    quiver(xarrow,yarrow,uarrow,varrow,'k','LineWidth',3.5)  % Plot shear vector
    hold off
    pbaspect([1 1 1])
    frameid=sprintf('%02s_RH_%03d_to_%03d_%02s_%02d_%02s_%02d.png',dset,lev(lev_rh(end)),lev(lev_rh(1)),tcplotname,t,shearplotname,s);
    set(gcf,'inverthardcopy','off')
    print(h,'-dpng',frameid);
    clf

  % FINALLY, SAME AS ABOVE, BUT FOR THE AZIMUTHAL ANOMALY OF RELATIVE HUMIDITY
    set(groot,'DefaultFigureColormap',flipud(parula))
    h=figure;
    imagesc(-5:xgrid:5,-5:ygrid:5,v5_new(:,:,t,s).',[min(min(v5_new(:,:,t,s))) max(max(v5_new(:,:,t,s)))]);
    hold on
    for c=1:5   % Plot circles to better visualize radial distances from center.
      viscircles([0 0],c,'Color','k','LineStyle','--','EnhanceVisibility',false,'LineWidth',1);
    end
    plot([0 0],[-5 5],'--k','LineWidth',1.0);
    plot([-5 5],[0 0],'--k','LineWidth',1.0);  % W in 0.5 Pa/s intervals (ascent only); theta anomalies in 0.2 K intervals
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-3 -3],'LineColor','k','LineStyle','-','LineWidth',2.5)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-2.5 -2.5],'LineColor','k','LineStyle','-','LineWidth',2.2)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-2 -2],'LineColor','k','LineStyle','-','LineWidth',1.9)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-1.5 -1.5],'LineColor','k','LineStyle','-','LineWidth',1.6)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-1 -1],'LineColor','k','LineStyle','-','LineWidth',1.3)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-0.5 -0.5],'LineColor','k','LineStyle','-','LineWidth',1.0)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[-1 -1],'LineColor','b','LineStyle','-','LineWidth',2.5)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[-0.8 -0.8],'LineColor','b','LineStyle','-','LineWidth',2.2)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[-0.6 -0.6],'LineColor','b','LineStyle','-','LineWidth',1.9)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[-0.4 -0.4],'LineColor','b','LineStyle','-','LineWidth',1.6)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[-0.2 -0.2],'LineColor','b','LineStyle','-','LineWidth',1.3)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[0.2 0.2],'LineColor','r','LineStyle','-','LineWidth',1.3)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[0.4 0.4],'LineColor','r','LineStyle','-','LineWidth',1.6)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[0.6 0.6],'LineColor','r','LineStyle','-','LineWidth',1.9)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[0.8 0.8],'LineColor','r','LineStyle','-','LineWidth',2.2)
    contour(lonbox,latbox,v3_new(:,:,t,s).',[1 1],'LineColor','r','LineStyle','-','LineWidth',2.5)
    xlim([-5 5])
    ylim([-5 5])
    set(gca,'XTick',[-5:1:5])
    set(gca,'YTick',[-5:1:5])
    set(gca,'YDir','normal')
    set(gca,'FontSize',14)
    cbh=colorbar;
    caxis([min(min(v5_new(:,:,t,s))) max(max(v5_new(:,:,t,s)))])   % May adjust colorbar limits for readability
    xlabel('Across-Shear Degrees from Center')
    ylabel('Along-Shear Degrees from Center')
    str1=strcat(dset," ",lev(lev_rh(end)),"-",lev(lev_rh(1))," hPa Relative Humidity Anomaly (%)");  % All necessary info will be included in plot title, and n$
    str2=strcat(tcbintext(t),", ",shearbintext(s));
    title([str1,str2])
    xarrow=0; yarrow=0; uarrow=0; varrow=5;
    quiver(xarrow,yarrow,uarrow,varrow,'k','LineWidth',3.5)  % Plot shear vector
    hold off
    pbaspect([1 1 1])
    frameid=sprintf('%02s_RHAnom_%03d_to_%03d_%02s_%02d_%02s_%02d.png',dset,lev(lev_rh(end)),lev(lev_rh(1)),tcplotname,t,shearplotname,s);
    set(gcf,'inverthardcopy','off')
    print(h,'-dpng',frameid);
    clf
  end
end
