% FUNCTION RainRate - Edited by Jake Carstens 10/12/2023 for manuscript on tropical cyclone asymmetry in reanalyses
%                     Generalized version of script used to generate Figures 5 and 7 of manuscript.
% PURPOSE: Compute and plot convective, large-scale, and total rain rate in shear-relative space, with contours of
%          vertical velocity, convergence, and CAPE overlaid. The vertical levels for the latter 3 variables may be
%          adjusted in the section where the variables are initially loaded.
% NOTES: 1. Storm-centered data snapshots, TempestExtremes-derived TC tracks, and wind shear files will be provided
%           in same directory within repository, so filepaths should not need to be altered.
%        2. The user has the choice to use either .mat or .nc format for the wind shear and storm-centered data. TempestExtremes
%           tracks are provided as .txt files.
% PROCEDURE:
% 1. Choose reanalysis dataset (currently ERA5 or CFSR), loading in its TC tracks, associated wind shear,
%    and information about the grid.
% 2. Pre-allocate matrices for convective (v1), large-scale (v2), and total (v3) rain rate, as well as vertical velocity
%    (v4), convergence (v5), and CAPE (v6), binned by both TC behavior and wind shear magnitude.
% 3. Loop through all TC snapshots, rotate data to shear vector, and place into TC/shear bin composite.
% 4. Test radial profiles for statistically significant differences between quadrant pairs, and plot each composite
%    radial profile and spatial map. Rain rate will be shaded, while other dynamic/thermodynamic variables will be contoured.

% FIRST, CHOOSE THE DATASET. THIS INCLUDES INFORMATION ABOUT THE GRID SPACING, USEFUL FOR DRAWING A 10-DEGREE BOX
% AROUND THE TC (x10d/y10d), AS WELL AS RE-CENTERING (xtotal) THE DOMAIN ABOUT THE TC.
dset="CFSR";
if (dset == "CFSR")
  x10d=10; y10d=10; xtotal=720; xgrid=0.50; ygrid=0.50; dset_c='CFSR';
elseif (dset == "ERA5")
  x10d=20; y10d=20; xtotal=1440; xgrid=0.25; ygrid=0.25; dset_c='ERA5';
end
lev=[100:25:250 300:50:750 775:25:1000];   % Vertical levels, which are the same in both reanalyses here. 27 total.
w_lev=21; conv_lev=25; % CHOOSE VERTICAL LEVELS TO DISPLAY VERTICAL VELOCITY AND CONVERGENCE DATA AT!
                       % CURRENTLY SET TO 850 hPa (w) AND 950 hPa (conv).

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
dist=[0:25:500];                                                        % Radii for quadrant profiles will be rounded to nearest 25 km (~0.25 deg).
count=zeros(length(tcbins)-1,length(shearbins)-1);                      % Accumulate all snapshots in a particular TC/shear bin.
v1_comp=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);  % Draw a 10-deg by 10-deg box for the TC composite.
v2_comp=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);
v3_comp=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);
v4_comp=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);  % I will choose the pressure levels for v4/v5 below when loading data.
v5_comp=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);
v6_comp=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);
v1_quad=NaN(4,length(dist),length(tcbins)-1,length(shearbins)-1,2500);  % Quadrant radial profiles. Order: DR --> DL --> UL --> UR
v2_quad=NaN(4,length(dist),length(tcbins)-1,length(shearbins)-1,2500);  % These will be kept discrete, in order to perform statistical significance
v3_quad=NaN(4,length(dist),length(tcbins)-1,length(shearbins)-1,2500);  % testing using the full distributions. (Only need this for rain rate!)

% IF CHOOSING THE .mat OPTION, QUEUE UP THE FILES HERE
CRR=matfile(['StormCenteredData/' dset_c '_CRR.mat']);
LSRR=matfile(['StormCenteredData/' dset_c '_LSRR.mat']);
W=matfile(['StormCenteredData/' dset_c '_W.mat']);
U=matfile(['StormCenteredData/' dset_c '_U.mat']);
V=matfile(['StormCenteredData/' dset_c '_V.mat']);
MLCAPE=matfile(['StormCenteredData/' dset_c '_MLCAPE.mat']);  % Using mixed-layer here, also have MU and model-output CAPE.

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
  u=U.u(:,:,:,p); % U wind
  v=V.v(:,:,:,p); % V wind
  w=W.w(:,:,:,p); % W wind
  crr=CRR.crr(:,:,p); % Convective parameterization rain rate
  lsrr=LSRR.lsrr(:,:,p); % Large-scale cloud scheme rain rate
  mlcape=MLCAPE.mlcape(:,:,p); % 0-1 km mixed layer CAPE (again, can be changed)

% .nc SECTION
  start=[1 1 1 p]; interval=[length([-10:xgrid:10]) length([-10:ygrid:10]) length(lev) 1];
  u=ncread(['StormCenteredData/' dset_c '_U.nc'],'u',start,interval);
  v=ncread(['StormCenteredData/' dset_c '_V.nc'],'v',start,interval);
  w=ncread(['StormCenteredData/' dset_c '_W.nc'],'w',start,interval);
  crr=ncread(['StormCenteredData/' dset_c '_CRR.nc'],'crr',start,interval);
  lsrr=ncread(['StormCenteredData/' dset_c '_LSRR.nc'],'lsrr',start,interval);
  mlcape=ncread(['StormCenteredData/' dset_c '_MLCAPE.nc'],'mlcape',start,interval);

  v1=crr; v2=lsrr; v3=v1+v2;   % ORDER WE WANT: CV --> LS --> TOTAL
  v4=w(:,:,w_lev);    % Vertical velocity (CHOOSE VERTICAL LEVEL AT TOP OF FUNCTION! Higher number --> higher pressure)
  v5=zeros(4*x10d+1,4*y10d+1); lat_box=lat(p)-10:ygrid:lat(p)+10; lon_box=lon(p)-10:xgrid:lon(p)+10;
  for a=2:4*x10d    % COMPUTE CONVERGENCE FROM U/V WINDS AT LEVEL SPECIFIED AT TOP
    for b=2:4*y10d
      dx=lldistkm([lat_box(b) lon_box(a-1)],[lat_box(b) lon_box(a+1)]).*1000.0; % Zonal distance between grid points
      dy=lldistkm([lat_box(b-1) lon_box(a)],[lat_box(b+1) lon_box(a)]).*1000.0; % Meridional distance (in m)
      du=u(a+1,b,conv_lev)-u(a-1,b,conv_lev);
      dv=v(a,b+1,conv_lev)-v(a,b-1,conv_lev);
      v5(a,b)=(du./dx)+(dv./dy);    % CONVERGENCE (or really, divergence. But we'll plot the negative contours later)
    end
  end
  v6=mlcape;

% USE A ROTATION MATRIX TO REORIENT EACH VARIABLE RELATIVE TO ITS SHEAR VECTOR.
% THIS IS DESIGNED TO HAVE THE SHEAR VECTOR POINTING UPWARDS ON A PLOT. (downshear --> positive y)
  dwnshr=atan2d(shear_v(p),shear_u(p));
  v1R=imrotate(v1,dwnshr.*(-1)+90,'nearest','crop');  % By extracting a larger box initially, we don't accidentally
  v2R=imrotate(v2,dwnshr.*(-1)+90,'nearest','crop');  % crop off data that we'll need when considering a 500 km area around the TC.
  v3R=imrotate(v3,dwnshr.*(-1)+90,'nearest','crop');
  v4R=imrotate(v4,dwnshr.*(-1)+90,'nearest','crop');
  v5R=imrotate(v5,dwnshr.*(-1)+90,'nearest','crop');
  v6R=imrotate(v6,dwnshr.*(-1)+90,'nearest','crop');
  [x,y]=meshgrid(lon(p)-(2*x10d*xgrid):xgrid:lon(p)+(2*x10d*xgrid),lat(p)-(2*y10d*ygrid):ygrid:lat(p)+(2*y10d*ygrid));
  latR=imrotate(y,dwnshr.*(-1)+90,'nearest','crop');  % Rotate the lats/lons as well, to make sure our distances are correct.
  lonR=imrotate(x,dwnshr.*(-1)+90,'nearest','crop');
  clear x y
  xc=floor(length(xint)/2)+1; yc=floor(length(yint)/2)+1;  % Extract the center of the new rotated square as the TC center to separate quadrants.

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
        count(t,s)=count(t,s)+1;        % Accumulate each snapshot falling into the bin, since we'll calculate a composite mean later.

% NOW THAT WE HAVE THE FULL BOX ACCOUNTED FOR, LET'S CONSIDER THE QUADRANT-SPECIFIC RADIAL PROFILES.
% Recall that I've only set this to consider the rain rate. One could simply add in CAPE, OMEGA, etc. using the same format.
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
          distcount=zeros(length(dist),1);
          for x=1:x10d+1
            for y=1:y10d+1     % Assign each point within the quadrant to a radius bin.
              distances=round(lldistkm([lat(p) lon(p)],[lat_quad(x,y) lon_quad(x,y)])/25)*25;  % Rounds to nearest 25 km
              v1_q(find(dist==distances))=v1_q(find(dist==distances))+v1R(xrange(x),yrange(y));
              v2_q(find(dist==distances))=v2_q(find(dist==distances))+v2R(xrange(x),yrange(y));
              v3_q(find(dist==distances))=v3_q(find(dist==distances))+v3R(xrange(x),yrange(y));
              distcount(find(dist==distances))=distcount(find(dist==distances))+1;
            end
          end
          distcount(distcount==0)=NaN;   % If any radii are unaccounted for, change the value to NaN.
          v1_quad(q,:,t,s,count(t,s))=v1_q./distcount;  % Take the mean of all points at similar radii within the snapshot
          v2_quad(q,:,t,s,count(t,s))=v2_q./distcount;  % but keep the snapshot itself discrete in order to conduct
          v3_quad(q,:,t,s,count(t,s))=v3_q./distcount;  % statistical significance testing later on.
        end
        break;
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
  end
end

% REWRITE CAPE AS THE ANOMALY FROM THE AZIMUTHAL MEAN
v6_azi=zeros(2*x10d+1,2*y10d+1,length(tcbins)-1,length(shearbins)-1);
distancesnew=zeros(2*x10d+1,2*y10d+1); xc=floor((2*x10d+1)/2)+1; yc=floor((2*y10d+1)/2)+1;
for a=1:2*x10d+1
  for b=1:2*y10d+1
    distancesnew(a,b)=sqrt((xc-a).^2+(yc-b).^2);
  end
end
for t=1:length(tcbins)-1
  for s=1:length(shearbins)-1
    v6_bin=reshape(v6_comp(:,:,t,s),2*x10d+1,2*y10d+1);
    for a=1:2*x10d+1
      for b=1:2*y10d+1
        distnew=sqrt((xc-a).^2+(yc-b).^2);
        v6_azi(a,b,t,s)=v6_bin(a,b)-mean(mean(v6_bin(find(distancesnew==distnew))));
      end
    end
  end
end

% PLOTTING SECTION - CHANGE PLOT STORAGE DIRECTORY ACCORDINGLY.
cd Plots
set(groot,'DefaultFigureColor','white')
precip_colors=[ 255 255 255  % RADAR-LIKE CUSTOM COLORMAP FOR SHADING RAIN RATE
 235 246 255
 214 226 255
 181 201 255
 142 178 255
 127 150 255
 114 133 248
  99 112 248
   0 158  30
  60 188  61
 185 249 110
 255 249  19
 255 163   9
 255  25  25
 204  61  61
 165  49  49
 237   0 237
 237 150 237]./255.0;

% PLOT EACH BIN OF SHEAR-RELATIVE VALUES, LOOKING DOWNSHEAR IN THE 10-DEGREE BOX
for t=1:length(tcbins)-1
  for s=1:length(shearbins)-1
    for r=1:length(dist)
%   THIS BLOCK SIMPLY TESTS EACH QUADRANT PAIR FOR STATISTICALLY SIGNIFICANT DIFFERENCES.
%   This outputs a corresponding p-value considering all points in the sample at a given
%   TC/shear/radius bin. I have explored avenues to indicate statistical significance
%   on the radial profile plot below, but have opted to not include that here.
%   This section only considers the rain rate variables. Others may be included by the user
%   with the same formatting (CAPE, convergence, omega).
      [n,v1_drdl]=ttest(v1_quad(1,r,s,1:count(t,s)),v1_quad(2,r,s,1:count(t,s)));
      [n,v2_drdl]=ttest(v2_quad(1,r,s,1:count(t,s)),v2_quad(2,r,s,1:count(t,s)));
      [n,v3_drdl]=ttest(v3_quad(1,r,s,1:count(t,s)),v3_quad(2,r,s,1:count(t,s)));
      [n,v1_drul]=ttest(v1_quad(1,r,s,1:count(t,s)),v1_quad(3,r,s,1:count(t,s)));
      [n,v2_drul]=ttest(v2_quad(1,r,s,1:count(t,s)),v2_quad(3,r,s,1:count(t,s)));
      [n,v3_drul]=ttest(v3_quad(1,r,s,1:count(t,s)),v3_quad(3,r,s,1:count(t,s)));
      [n,v1_drur]=ttest(v1_quad(1,r,s,1:count(t,s)),v1_quad(4,r,s,1:count(t,s)));
      [n,v2_drur]=ttest(v2_quad(1,r,s,1:count(t,s)),v2_quad(4,r,s,1:count(t,s)));
      [n,v3_drur]=ttest(v3_quad(1,r,s,1:count(t,s)),v3_quad(4,r,s,1:count(t,s)));
      [n,v1_dlul]=ttest(v1_quad(2,r,s,1:count(t,s)),v1_quad(3,r,s,1:count(t,s)));
      [n,v2_dlul]=ttest(v2_quad(2,r,s,1:count(t,s)),v2_quad(3,r,s,1:count(t,s)));
      [n,v3_dlul]=ttest(v3_quad(2,r,s,1:count(t,s)),v3_quad(3,r,s,1:count(t,s)));
      [n,v1_dlur]=ttest(v1_quad(2,r,s,1:count(t,s)),v1_quad(4,r,s,1:count(t,s)));
      [n,v2_dlur]=ttest(v2_quad(2,r,s,1:count(t,s)),v2_quad(4,r,s,1:count(t,s)));
      [n,v3_dlur]=ttest(v3_quad(2,r,s,1:count(t,s)),v3_quad(4,r,s,1:count(t,s)));
      [n,v1_urul]=ttest(v1_quad(3,r,s,1:count(t,s)),v1_quad(4,r,s,1:count(t,s)));
      [n,v2_urul]=ttest(v2_quad(3,r,s,1:count(t,s)),v2_quad(4,r,s,1:count(t,s)));
      [n,v3_urul]=ttest(v3_quad(3,r,s,1:count(t,s)),v3_quad(4,r,s,1:count(t,s)));
    end

%   FIRST, PLOT THE QUADRANT RADIAL PROFILES. HERE, I OVERLAY CONVECTIVE/LARGE-SCALE RAIN ON ONE PLOT,
%   AND ONLY SHOW THE TOTAL RAIN RATE ON ANOTHER FOR READABILITY.

    c_dr=[254 97 0]./255; c_dl=[255 176 0]./255; c_ul=[100 143 255]./255; c_ur=[220 38 127]./255;
    s_v1='^'; s_v2='s'; s_v3='o'; sz=20;  % Simply set the symbols/colors here for different variables/quadrants.

    h=figure;  % nanmean will make sure that any blank data points will not be considered.
    plot(dist,nanmean(v1_quad(1,:,t,s,1:count(s)),5),'Color',c_dr,'Marker',s_v1,'LineWidth',1.8);
    hold on    % Take the mean in the 5th dimension, or across the sample at a given radius/shear/TC characteristic.
    plot(dist,nanmean(v1_quad(2,:,t,s,1:count(s)),5),'Color',c_dl,'Marker',s_v1,'LineWidth',1.8);
    plot(dist,nanmean(v1_quad(3,:,t,s,1:count(s)),5),'Color',c_ul,'Marker',s_v1,'LineWidth',1.8);
    plot(dist,nanmean(v1_quad(4,:,t,s,1:count(s)),5),'Color',c_ur,'Marker',s_v1,'LineWidth',1.8);
    plot(dist,nanmean(v2_quad(1,:,t,s,1:count(s)),5),'Color',c_dr,'Marker',s_v2,'LineWidth',1.8);
    plot(dist,nanmean(v2_quad(2,:,t,s,1:count(s)),5),'Color',c_dl,'Marker',s_v2,'LineWidth',1.8);
    plot(dist,nanmean(v2_quad(3,:,t,s,1:count(s)),5),'Color',c_ul,'Marker',s_v2,'LineWidth',1.8);
    plot(dist,nanmean(v2_quad(4,:,t,s,1:count(s)),5),'Color',c_ur,'Marker',s_v2,'LineWidth',1.8);
    hold off
    grid on
    xlim([0 max(dist)])     % As it stands, limit the plot to a 0-500 km radial range
%    ylim([0 max(max(max(max(max(v2_quad)))))+1])  % May be adjusted for readability (large-scale rain usually larger)
    set(gca,'FontSize',14)                      % For now, I'll have the y-axis scaling with the data in each bin.
    xlabel('Distance (km)')
    ylabel('Rain Rate (mm/hr)')
    str1=strcat(dset," Decomposed Rain Rate");
    str2=strcat(tcbintext(t),", ",shearbintext(s));    % Write out details of compositing in plot title.
    title([str1,str2])
    leg=legend('DR - CV','DL - CV','UL - CV','UR - CV','DR - LS','DL - LS','UL - LS','UR - LS','location','eastoutside');
    leg.FontSize=11;
    hold off
    frameid=sprintf('%02s_QuadrantRainDecomposed_%02s_%02d_%02s_%02d.png',dset,tcplotname,t,shearplotname,s);   % Generates file name.
    set(gcf,'inverthardcopy','off')
    print(h,'-dpng',frameid);
    clf

    h=figure;  % nanmean will make sure that any blank data points will not be considered.
    plot(dist,nanmean(v3_quad(1,:,t,s,1:count(s)),5),'Color',c_dr,'Marker',s_v3,'LineWidth',1.8);
    hold on    % Take the mean in the 5th dimension, or across the sample at a given radius/shear/TC characteristic.
    plot(dist,nanmean(v3_quad(2,:,t,s,1:count(s)),5),'Color',c_dl,'Marker',s_v3,'LineWidth',1.8);
    plot(dist,nanmean(v3_quad(3,:,t,s,1:count(s)),5),'Color',c_ul,'Marker',s_v3,'LineWidth',1.8);
    plot(dist,nanmean(v3_quad(4,:,t,s,1:count(s)),5),'Color',c_ur,'Marker',s_v3,'LineWidth',1.8);
    hold off
    grid on
    xlim([0 max(dist)])     % As it stands, limit the plot to a 0-500 km radial range
%    ylim([0 max(max(max(max(max(v3_quad)))))+1])  % May be adjusted for readability (large-scale rain usually larger)
    set(gca,'FontSize',14)                      % For now, I'll have the y-axis scaling with the data in each bin.
    xlabel('Distance (km)')
    ylabel('Rain Rate (mm/hr)')
    str1=strcat(dset," Total Rain Rate");
    str2=strcat(tcbintext(t),", ",shearbintext(s));    % Write out details of compositing in plot title.
    title([str1,str2])
    leg=legend('DR','DL','UL','UR','location','eastoutside');  % May adjust legend location to optimize figure size/data visibility.
    leg.FontSize=11;
    hold off
    frameid=sprintf('%02s_QuadrantRainTotal_%02s_%02d_%02s_%02d.png',dset,tcplotname,t,shearplotname,s);   % Generates file name.
    set(gcf,'inverthardcopy','off')
    print(h,'-dpng',frameid);
    clf

%   NEXT, PLOT A SPATIAL MAP OF TOTAL RAIN RATE WITH THE SHEAR VECTOR POINTING UPWARDS.
    gray=[0.5,0.5,0.5]; lonbox=-5:xgrid:5; latbox=-5:ygrid:5;
    set(groot,'DefaultFigureColormap',precip_colors)  % Sets colormap to custom version input above
    h=figure;
    imagesc(lonbox,latbox,v3_comp(:,:,t,s).',[0 max(max(max(max(v3_comp))))]);
    hold on
    for c=1:5  % Plot circles to show radii in 1-deg intervals
      viscircles([0 0],c,'Color','k','LineStyle','--','EnhanceVisibility',false,'LineWidth',1);
    end
    plot([0 0],[-5 5],'--k','LineWidth',1.0);
    plot([-5 5],[0 0],'--k','LineWidth',1.0);
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-0.25 -0.25],'LineColor',gray,'LineStyle','-','LineWidth',1.3)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-0.5 -0.5],'LineColor',gray,'LineStyle','-','LineWidth',1.6)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-1 -1],'LineColor',gray,'LineStyle','-','LineWidth',1.9)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-1.5 -1.5],'LineColor',gray,'LineStyle','-','LineWidth',1.2)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-2 -2],'LineColor',gray,'LineStyle','-','LineWidth',2.5)
    xlim([-5 5])    % Contours above show vertical velocity, with thicker contours --> stronger ascent
    ylim([-5 5])
    set(gca,'XTick',[-5:1:5])
    set(gca,'YTick',[-5:1:5])
    set(gca,'YDir','normal')
    set(gca,'FontSize',14)
    cbh=colorbar;
    caxis([0 max(max(max(max(v3_comp))))]) % May adjust colorbar limits for readability
    xlabel('Across-Shear Degrees from Center')
    ylabel('Along-Shear Degrees from Center')
    str1=strcat(dset," Total Rain Rate (mm/hr)");
    str2=strcat(tcbintext(t),", ",shearbintext(s));
    title([str1,str2])
    xarrow=0; yarrow=0; uarrow=0; varrow=5;
    quiver(xarrow,yarrow,uarrow,varrow,'k','LineWidth',3.5)  % Plot shear vector
    hold off
    pbaspect([1 1 1])
    frameid=sprintf('%02s_TotalRain_%02s_%02d_%02s_%02d.png',dset,tcplotname,t,shearplotname,s);
    set(gcf,'inverthardcopy','off')
    print(h,'-dpng',frameid);
    clf

% NEXT UP, THE CONVECTIVE RAIN FOLLOWING A SIMILAR PROCEDURE. WE'LL ALSO OVERLAY CONVERGENCE AND CAPE
% CONTOURS HERE, USING THE AZIMUTHAL CAPE ANOMALY (v6_new) AS COMPUTED ABOVE.
    set(groot,'DefaultFigureColormap',precip_colors)
    h=figure;
    imagesc(lonbox,latbox,v1_comp(:,:,t,s).',[0 max(max(max(max(v1_comp))))]);
    hold on
    for c=1:5
      viscircles([0 0],c,'Color','k','LineStyle','--','EnhanceVisibility',false,'LineWidth',1);
    end
    plot([0 0],[-5 5],'--k','LineWidth',1.0);
    plot([-5 5],[0 0],'--k','LineWidth',1.0);
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-0.25 -0.25],'LineColor',gray,'LineStyle','-','LineWidth',1.3)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-0.5 -0.5],'LineColor',gray,'LineStyle','-','LineWidth',1.6)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-1 -1],'LineColor',gray,'LineStyle','-','LineWidth',1.9)  % Omega
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-1.5 -1.5],'LineColor',gray,'LineStyle','-','LineWidth',2.2)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-2 -2],'LineColor',gray,'LineStyle','-','LineWidth',2.5)
    contour(lonbox,latbox,v5_comp(:,:,t,s).',[-0.00003 -0.00003],'LineColor','b','LineStyle','-','LineWidth',1.3) % Convergence
    contour(lonbox,latbox,v5_comp(:,:,t,s).',[-0.00006 -0.00006],'LineColor','b','LineStyle','-','LineWidth',1.6)
    contour(lonbox,latbox,v5_comp(:,:,t,s).',[-0.00009 -0.00009],'LineColor','b','LineStyle','-','LineWidth',1.9)
    contour(lonbox,latbox,v5_comp(:,:,t,s).',[-0.00012 -0.00012],'LineColor','b','LineStyle','-','LineWidth',2.2)
    contour(lonbox,latbox,v6_azi(:,:,t,s).',[50 50],'LineColor','r','LineStyle','-','LineWidth',1.3)  % CAPE azimuthal anomaly
    contour(lonbox,latbox,v6_azi(:,:,t,s).',[100 100],'LineColor','r','LineStyle','-','LineWidth',1.6)
    contour(lonbox,latbox,v6_azi(:,:,t,s).',[150 150],'LineColor','r','LineStyle','-','LineWidth',1.9)
    contour(lonbox,latbox,v6_azi(:,:,t,s).',[200 200],'LineColor','r','LineStyle','-','LineWidth',2.2)
    contour(lonbox,latbox,v6_azi(:,:,t,s).',[250 250],'LineColor','r','LineStyle','-','LineWidth',2.5)
    xlim([-5 5])
    ylim([-5 5])
    set(gca,'XTick',[-5:1:5])
    set(gca,'YTick',[-5:1:5])
    set(gca,'YDir','normal')
    set(gca,'FontSize',14)
    cbh=colorbar;
    caxis([0 max(max(max(max(v1_comp))))])
    xlabel('Across-Shear Degrees from Center')
    ylabel('Along-Shear Degrees from Center')
    str1=strcat(dset," Convective Rain (mm/hr)");
    str2=strcat(tcbintext(t),", ",shearbintext(s));
    title([str1,str2])
    xarrow=0; yarrow=0; uarrow=0; varrow=5;
    quiver(xarrow,yarrow,uarrow,varrow,'k','LineWidth',3.5)
    hold off
    pbaspect([1 1 1])
    frameid=sprintf('%02s_ConvectiveRain_%02s_%02d_%02s_%02d.png',dset,tcplotname,t,shearplotname,s);
    set(gcf,'inverthardcopy','off')
    print(h,'-dpng',frameid);
    clf

% AND FINALLY, THE LARGE-SCALE RAIN RATE
    set(groot,'DefaultFigureColormap',precip_colors)
    h=figure;
    imagesc(lonbox,latbox,v2_comp(:,:,t,s).',[0 max(max(max(max(v2_comp))))]);
    hold on
    for c=1:5
      viscircles([0 0],c,'Color','k','LineStyle','--','EnhanceVisibility',false,'LineWidth',1);
    end
    plot([0 0],[-5 5],'--k','LineWidth',1.0);
    plot([-5 5],[0 0],'--k','LineWidth',1.0);
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-0.25 -0.25],'LineColor',gray,'LineStyle','-','LineWidth',1.3)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-0.5 -0.5],'LineColor',gray,'LineStyle','-','LineWidth',1.6)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-1 -1],'LineColor',gray,'LineStyle','-','LineWidth',1.9)  % Omega
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-1.5 -1.5],'LineColor',gray,'LineStyle','-','LineWidth',2.2)
    contour(lonbox,latbox,v4_comp(:,:,t,s).',[-2 -2],'LineColor',gray,'LineStyle','-','LineWidth',2.5)
    contour(lonbox,latbox,v5_comp(:,:,t,s).',[-0.00003 -0.00003],'LineColor','b','LineStyle','-','LineWidth',1.3) % Convergence
    contour(lonbox,latbox,v5_comp(:,:,t,s).',[-0.00006 -0.00006],'LineColor','b','LineStyle','-','LineWidth',1.6)
    contour(lonbox,latbox,v5_comp(:,:,t,s).',[-0.00009 -0.00009],'LineColor','b','LineStyle','-','LineWidth',1.9)
    contour(lonbox,latbox,v5_comp(:,:,t,s).',[-0.00012 -0.00012],'LineColor','b','LineStyle','-','LineWidth',2.2)
    contour(lonbox,latbox,v6_azi(:,:,t,s).',[50 50],'LineColor','r','LineStyle','-','LineWidth',1.3)  % CAPE azimuthal anomaly
    contour(lonbox,latbox,v6_azi(:,:,t,s).',[100 100],'LineColor','r','LineStyle','-','LineWidth',1.6)
    contour(lonbox,latbox,v6_azi(:,:,t,s).',[150 150],'LineColor','r','LineStyle','-','LineWidth',1.9)
    contour(lonbox,latbox,v6_azi(:,:,t,s).',[200 200],'LineColor','r','LineStyle','-','LineWidth',2.2)
    contour(lonbox,latbox,v6_azi(:,:,t,s).',[250 250],'LineColor','r','LineStyle','-','LineWidth',2.5)
    xlim([-5 5])
    ylim([-5 5])
    set(gca,'XTick',[-5:1:5])
    set(gca,'YTick',[-5:1:5])
    set(gca,'YDir','normal')
    set(gca,'FontSize',14)
    cbh=colorbar;
    caxis([0 max(max(max(max(v2_comp))))])
    xlabel('Across-Shear Degrees from Center')
    ylabel('Along-Shear Degrees from Center')
    str1=strcat(dset," Large-Scale Rain (mm/hr)");
    str2=strcat(tcbintext(t),", ",shearbintext(s));
    title([str1,str2])
    xarrow=0; yarrow=0; uarrow=0; varrow=5;
    quiver(xarrow,yarrow,uarrow,varrow,'k','LineWidth',3.5)
    hold off
    pbaspect([1 1 1])
    frameid=sprintf('%02s_LargeScaleRain_%02s_%02d_%02s_%02d.png',dset,tcplotname,t,shearplotname,s);
    set(gcf,'inverthardcopy','off')
    print(h,'-dpng',frameid);
    clf
  end
end
