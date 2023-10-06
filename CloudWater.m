% FUNCTION CloudWater - Edited by Jake Carstens 10/6/2023 for manuscript on tropical cyclone asymmetry in reanalyses
%                       Generalized version of script used to generate Figure 8 of manuscript.
% PURPOSE: Compute and plot cloud liquid and ice water content in reanalyses, vertically-resolved across 4 shear-relative
%          quadrants. Overlay various kinematic and thermodynamic properties on top of this as contours, including
%          vertical motion, convergence, and the freezing/melting level.
% NOTES: 1. Storm-centered data snapshots, TempestExtremes-derived TC tracks, and wind shear files will be provided
%           in same directory within repository, so filepaths should not need to be altered.
%        2. The user has the choice to use either .mat or .nc format for the wind shear and storm-centered data. TempestExtremes
%           tracks are provided as .txt files.
%        3. CFSR DOES NOT PARTITION CLOUD LIQUID AND ICE! For CFSR, the "clw" variable represents both liquid and ice combined.
% PROCEDURE:
% 1. Load in TC tracks, associated wind shear, and information about the grid.
% 2. Pre-allocate matrices for cloud liquid (v1), cloud ice (v2), vertical velocity (v3), divergence (v4), radial wind
%    (v5), and temperature (v6), binned by both TC behavior and wind shear magnitude. r/z profiles for each quadrant.
% 3. Loop through all TC snapshots, rotate data to shear vector, and place into TC/shear bin composite.
% 4. Plot each composite radius/pressure profile of all variables overlaid, for each quadrant separately.

% FIRST, STATE WHAT REANALYSIS/MODEL GRID WE'RE USING TO DEFINE THE 10-DEGREE BOXES
dset="CFSR";
if (dset == "CFSR")
  x10d=10; y10d=10; xtotal=720; xgrid=0.50; ygrid=0.50; dset_c='CFSR'; dset_n='CFSR';
elseif (dset == "ERA5")
  x10d=20; y10d=20; xtotal=1440; xgrid=0.25; ygrid=0.25; dset_c='ERA5'; dset="ERA5";
end
lev=[100:25:250 300:50:750 775:25:1000];   % 27 total vertical levels.

% LOAD IN TC TRACKS FROM TEMPESTEXTREMES. LOAD IN STORM ID, TIME, POSITION, WIND SPEED, AND MINIMUM PRESSURE.
A=readmatrix(['TCTracks/trajectories_' dset_c '.txt']);
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
% YOU MAY ALSO CHANGE THE BOUNDARIES FOR THE BINS HERE!
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
% v1: CLOUD LIQUID  v2: CLOUD ICE  v3: OMEGA  v4: DIVERGENCE  v5: RADIAL WIND  v6: TEMPERATURE
dist=[0:25:500]; count=zeros(length(lev),length(tcbins)-1,length(shearbins)-1);  % Radii will be rounded to nearest 25 km (~0.25 deg).
v1_quad=NaN(4,length(dist),length(lev),length(tcbins)-1,length(shearbins)-1);    % Accumulate count falling into each TC/shear bin throughout.
if (dset == "ERA5")
  v2_quad=NaN(4,length(dist),length(lev),length(tcbins)-1,length(shearbins)-1);  % Quadrant radial profiles. Order: DR --> DL --> UL --> UR
end
v3_quad=NaN(4,length(dist),length(lev),length(tcbins)-1,length(shearbins)-1);
v4_quad=NaN(4,length(dist),length(lev),length(tcbins)-1,length(shearbins)-1);
v5_quad=NaN(4,length(dist),length(lev),length(tcbins)-1,length(shearbins)-1);
v6_quad=NaN(4,length(dist),length(lev),length(tcbins)-1,length(shearbins)-1);

% IF WORKING WITH .mat DATA, QUEUE UP THE FILES FOR ACCESS HERE. IF CHOOSING .nc, COMMENT THESE LINES OUT.
CLW=matfile(['StormCenteredData/' dset_c '_CLW.mat']);
if (dset == "ERA5")
  CIW=matfile(['StormCenteredData/' dset_c '_CIW.mat']);
end
W=matfile(['StormCenteredData/' dset_c '_W.mat']);
U=matfile(['StormCenteredData/' dset_c '_U.mat']);
V=matfile(['StormCenteredData/' dset_c '_V.mat']);
T=matfile(['StormCenteredData/' dset_c '_T.mat']);

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

% LOAD DATA HERE! BE SURE TO COMMENT OUT WHATEVER DATA FORMAT YOU ARE NOT USING BETWEEN .nc AND .mat!
% .mat SECTION
  clw=CLW.clw(:,:,:,p); w=W.w(:,:,:,p); u=U.u(:,:,:,p); v=V.v(:,:,:,p); temp=T.temp(:,:,:,p);
  if (dset == "ERA5")
    ciw=CIW.ciw(:,:,:,p);
  end
% .nc SECTION
%  start=[1 1 1 p]; interval=[length([-10:xgrid:10]) length([-10:ygrid:10]) length(lev) 1];
%  clw=ncread(['StormCenteredData/' dset_c '_CLW.nc'],'clw',start,interval);
%  if (dset == "ERA5")
%    ciw=ncread(['StormCenteredData/' dset_c '_CIW.nc'],'ciw',start,interval);
%  w=ncread(['StormCenteredData/' dset_c '_W.nc'],'w',start,interval);
%  u=ncread(['StormCenteredData/' dset_c '_U.nc'],'u',start,interval);
%  v=ncread(['StormCenteredData/' dset_c '_V.nc'],'v',start,interval);
%  temp=ncread(['StormCenteredData/' dset_c '_T.nc'],'temp',start,interval);
  % Need u and v winds to calculate divergence and radial wind!

  latbox=[-10:ygrid:10]+lat(p); lonbox=[-10:xgrid:10]+lon(p);
  v4=zeros(length(xint),length(yint),length(lev));
  for z=1:length(lev)   % Need to have tanradwinds.m utility in the same folder as this function!
    v1(:,:,z)=clw(:,:,z).*1000.0;   % CLOUD LIQUID (convert to g/kg)
    if (dset == "ERA5")
      v2(:,:,z)=ciw(:,:,z).*1000.0; % CLOUD ICE (convert to g/kg)
    end
    v3(:,:,z)=w(:,:,z);             % VERTICAL VELOCITY (Pa/s)
    [utan,urad]=tanradwinds(u(:,:,z),v(:,:,z));
    v5(:,:,z)=urad;                 % RADIAL WIND (m/s)
    v6(:,:,z)=temp(:,:,z)-273.15;   % TEMPERATURE (convert to C for simplicity)
    for a=2:4*x10d
      for b=2:4*y10d
        dx=lldistkm([lat_box(b) lon_box(a-1)],[lat_box(b) lon_box(a+1)]).*1000.0;
        dy=lldistkm([lat_box(b-1) lon_box(a)],[lat_box(b+1) lon_box(a)]).*1000.0;
        du=v_4(xint(a+1),yint(b),z,hr)-v_4(xint(a-1),yint(b),z,hr);
        dv=v_5(xint(a),yint(b+1),z,hr)-v_5(xint(a),yint(b-1),z,hr);
        v4(a,b,z)=(du./dx)+(dv./dy);       % DIVERGENCE
      end
    end
  end

% USE A ROTATION MATRIX TO REORIENT EACH VARIABLE RELATIVE TO ITS SHEAR VECTOR.
% THIS IS DESIGNED TO HAVE THE SHEAR VECTOR POINTING UPWARDS ON A PLOT. (downshear --> positive y)
% BY EXTRACTING A LARGER BOX AT FIRST, WE AVOID AN ISSUE OF CROPPING OUT DATA WE'D ACTUALLY NEED.
  for z=1:length(lev)
    dwnshr=atan2d(shear_v(p),shear_u(p));
    v1R=imrotate3(v1,dwnshr.*(-1)+90,[0 0 1],'nearest','crop');  % This rotates about the z-axis
    if (dset == "ERA5")
      v2R=imrotate3(v2,dwnshr.*(-1)+90,[0 0 1],'nearest','crop');
    end
    v3R=imrotate3(v3,dwnshr.*(-1)+90,[0 0 1],'nearest','crop');
    v4R=imrotate3(v4,dwnshr.*(-1)+90,[0 0 1],'nearest','crop');
    v5R=imrotate3(v5,dwnshr.*(-1)+90,[0 0 1],'nearest','crop');
    v6R=imrotate3(v6,dwnshr.*(-1)+90,[0 0 1],'nearest','crop');
    [x,y]=meshgrid(lon(p)-(2*x10d*xgrid):xgrid:lon(p)+(2*x10d*xgrid),lat(p)-(2*y10d*ygrid):ygrid:lat(p)+(2*y10d*ygrid));
    latR=imrotate(y,dwnshr.*(-1)+90,'nearest','crop');  % Rotate the lats/lons as well, to make sure our distances are correct.
    lonR=imrotate(x,dwnshr.*(-1)+90,'nearest','crop');
    clear x y
    xc=floor(length(xint)/2)+1; yc=floor(length(yint)/2)+1;  % Extract the center of the new rotated square as the TC center to separate quadrants.

% NOW, PLACE VARIABLES INTO BINS BASED ON THE SCHEMES PRESCRIBED INITIALLY
    for t=1:length(tcbins)-1
      for s=1:length(shearbins)-1
        if (shearvar(p) >= shearbins(s) && shearvar(p) < shearbins(s+1) && tcvar(p) >= tcbins(t) && tcvar(p) < tcbins(t+1))
          count(z,t,s)=count(z,t,s)+1;
          for q=1:4     % Order of quadrants: DR --> DL --> UL --> UR
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
            v1_q=zeros(length(dist),1); v3_q=zeros(length(dist),1);
            v4_q=zeros(length(dist),1); v5_q=zeros(length(dist),1); v6_q=zeros(length(dist),1);
            if (dset == "ERA5")
              v2_q=zeros(length(dist),1);
            end
            distcount=zeros(length(dist),1);
            for x=1:x10d+1
              for y=1:y10d+1     % Assign each point within the quadrant to a radius bin.
                distances=round(lldistkm([lat(p) lon(p)],[lat_quad(x,y) lon_quad(x,y)])/25)*25;  % Rounds to nearest 25 km
                v1_q(find(dist==distances))=v1_q(find(dist==distances))+v1R(xrange(x),yrange(y),z);
                if (dset == "ERA5")
                  v2_q(find(dist==distances))=v2_q(find(dist==distances))+v2R(xrange(x),yrange(y),z);
                end
                v3_q(find(dist==distances))=v3_q(find(dist==distances))+v3R(xrange(x),yrange(y),z);
                v4_q(find(dist==distances))=v4_q(find(dist==distances))+v4R(xrange(x),yrange(y),z);
                v5_q(find(dist==distances))=v5_q(find(dist==distances))+v5R(xrange(x),yrange(y),z);
                v6_q(find(dist==distances))=v6_q(find(dist==distances))+v6R(xrange(x),yrange(y),z);
                distcount(find(dist==distances))=distcount(find(dist==distances))+1;
              end
            end
            distcount(distcount==0)=NaN;   % If any radii are unaccounted for, change the value to NaN.
            v1_quad(q,:,z,t,s,count(z,t,s))=v1_q./distcount;    % Take the mean of all points at similar radii within the snapshot
            if (dset == "ERA5")
              v2_quad(q,:,z,t,s,count(z,t,s))=v2_q./distcount;
            end
            v3_quad(q,:,z,t,s,count(z,t,s))=v3_q./distcount;
            v4_quad(q,:,z,t,s,count(z,t,s))=v4_q./distcount;
            v5_quad(q,:,z,t,s,count(z,t,s))=v5_q./distcount;
            v6_quad(q,:,z,t,s,count(z,t,s))=v6_q./distcount;
          end
          break;
        end
      end
    end
  end
  toc
end

% PLOTTING SECTION - CHANGE PLOT STORAGE DIRECTORY ACCORDINGLY.
cd Plots   % Change to whatever directory you would like your plots stored in your workspace.
set(groot,'DefaultFigureColor','white')
set(groot,'DefaultFigureColormap',flipud(gray))  % Change to your preferred colormap (keeping it to MATLAB defaults here!)

% PLOT EACH BIN OF SHEAR-RELATIVE VALUES, LOOKING DOWNSHEAR IN THE 10-DEGREE BOX
for t=1:length(tcbins)-1
  for s=1:length(shearbins)-1
    for q=1:4
      if (q == 1)
        quadrant="Downshear Right Quadrant";
      elseif (q == 2)
        quadrant="Downshear Left Quadrant";
      elseif (q == 3)
        quadrant="Upshear Left Quadrant";
      elseif (q == 4)
        quadrant="Upshear Right Quadrant";
      end

      h=figure;
      contourf(dist,lev,nanmean(v1_quad(q,:,:,t,s,:),6).','LineColor','none','LevelStep',0.001)   % Shaded cloud liquid
      hold on
      if (dset == "ERA5")
        contour(dist,lev,nanmean(v2_quad(q,:,:,t,s,:),6).',[0.03:0.03:0.3],'LineColor','k','LineStyle','-','LineWidth',1.5)  % Cloud ice in black contours
      end
      contour(dist,lev,nanmean(v3_quad(q,:,:,t,s,:),6).',[-2:0.2:-0.2],'LineColor','r','LineStyle','-','LineWidth',1.5)    % Vertical velocity in red
      contour(dist,lev,nanmean(v3_quad(q,:,:,t,s,:),6).',[0.2:0.2:2],'LineColor','r','LineStyle','--','LineWidth',1.5)   % Descent in dashed lines
      contour(dist,lev,nanmean(v4_quad(q,:,:,t,s,:),6).',[0.00001:0.00002:0.00009],'LineColor','b','LineStyle','--','LineWidth',1.5)    % Divergence
      contour(dist,lev,nanmean(v4_quad(q,:,:,t,s,:),6).',[-0.00001:-0.00002:-0.00009],'LineColor','b','LineStyle','-','LineWidth',1.5)  % Convergence
      contour(dist,lev,nanmean(v6_quad(q,:,:,t,s,:),6).',[0 0],'LineColor','m','LineStyle','-','LineWidth',2.5)    % Freezing level as thick magenta contour
      quiver(dist,lev,nanmean(v5_quad(q,:,:,t,s,:),6).',nanmean(v3_quad(q,:,:,t,s,:),6).','b')      % Vectors of radial/vertical wind
      xlim([0 500])
      ylim([150 975])      % Can change visualization limits for aesthetics and closer examination of particular radii/vertical levels
      set(gca,'XTick',[0:50:500])
      set(gca,'YTick',[200:100:900])
      set(gca,'YDir','reverse')   % Simply sets pressure to decrease with height as in reality
      if (q == 2 || q == 3)
        set(gca,'XDir','reverse') % Flips x-axis to look more intuitively like we're looking left of shear
      else
        set(gca,'XDir','normal')
      end
      set(gca,'YScale','log')
      set(gca,'FontSize',14)
      colorbar;
      grid on
      xlabel('Radius (km)')
      ylabel('Pressure (hPa)')
      str1=strcat(dset," Cloud Water (g/kg)");
      str2=strcat(tcbintext(t),", ",shearbintext(s));
      str3=quadrant;
      title([str1,str2,str3])
      frameid=sprintf('%02s_CloudWater_%02s_%02d_%02s_%02d_%02d.png',dset,tcplotname,t,shearplotname,s,q);
      set(gcf,'inverthardcopy','off')
      print(h,'-dpng',frameid);
      clf
    end
  end
end
