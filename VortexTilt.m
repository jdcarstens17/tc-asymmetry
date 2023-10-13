% FUNCTION VortexTilt - Edited by Jake Carstens 5/10/2023 for manuscript on tropical cyclone asymmetry in reanalyses
%                       Generalized version of script used to generate Figure 4 of manuscript.
% PURPOSE: Calculate TC/shear-binned composite mean vortex tilt, relative to the surface pressure minimum, at regular
%          intervals above the surface based on a centroid of relative vorticity or geopotential height.
% NOTES: 1. Storm-centered data snapshots, TempestExtremes-derived TC tracks, and wind shear files will be provided
%           in same directory within repository, so filepaths should not need to be altered.
%        2. The user has the choice to use either .mat or .nc format for the wind shear and storm-centered data. TempestExtre$
%           tracks are provided as .txt files.
% PROCEDURE:
% 1. Choose reanalysis dataset (currently ERA5 or CFSR), loading in its TC tracks, associated wind shear,
%    and information about the grid.
% 2. Pre-allocate matrices for vortex center displacement in both x and y directions, binned by both TC behavior
%    (intensity, intensification rate, etc.) and wind shear magnitude.
% 3. Loop through all TC snapshots, rotate vorticity data to shear vector, calculate the zonal/meridional displacement
%    of the vorticity centroid at each vertical level relative to the surface center, and add to composite bin.
% 4. Plot each composite mean vortex tilt structure across all TC bins for a given shear bin, with the surface PMIN
%    at the center of the plot. For example, one plot will overlay the composite for TDs, weak TS, strong TS, and HU.

% FIRST, CHOOSE THE DATASET. THIS INCLUDES INFORMATION ABOUT THE GRID SPACING, USEFUL FOR DRAWING A 10-DEGREE BOX
% AROUND THE TC (x10d/y10d), AS WELL AS RE-CENTERING (xtotal) THE DOMAIN ABOUT THE TC.
dset="ERA5";
if (dset == "CFSR")
  x10d=10; y10d=10; xtotal=720; xgrid=0.50; ygrid=0.50; dset_c='CFSR';
elseif (dset == "ERA5")
  x10d=20; y10d=20; xtotal=1440; xgrid=0.25; ygrid=0.25; dset_c='ERA5';
end
lev=[100:25:250 300:50:750 775:25:1000];   % Vertical levels, which are the same in both reanalyses here. 27 total.
levs=[10:2:16 19 23];  % Plot centers from 400-900 hPa in 100 hPa intervals, with surface center (PMIN) at center of plot.

% THEN, CHOOSE THE METHOD BY WHICH YOU WANT TO ESTIMATE THE VORTEX CENTERS ALOFT.
% CURRENT OPTIONS ARE "VORT" (relative vorticity centroid) or "Z" (geopotential height centroid).
center_find_method="VORT";

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

% IF CHOOSING THE .mat OPTION, QUEUE UP THE FILES HERE
U=matfile(['StormCenteredData/' dset_c '_U.mat']);
V=matfile(['StormCenteredData/' dset_c '_V.mat']);
Z=matfile(['StormCenteredData/' dset_c '_Z.mat']);

% LOOP THROUGH SCHEMES TO BIN THE TCs. OPTIONS INCLUDE WIND SPEED, MINIMUM PRESSURE, AND WIND SPEED CHANGE IN THE NEXT 6 HOURS.
for tcbinscheme=1:3
  if (tcbinscheme == 1)  % Feel free to adjust bin boundaries as you wish.
    tcbins=[10 17 25 33 100]; tcbintext=["TD","Weak TS","Strong TS","Hurricane"];  % Winds are in m/s. Text will go in legend.
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

  count=zeros(length(levs),length(tcbins)-1,length(shearbins)-1);   % Keep track of snapshot counts to eventually calculate a composite mean.
  x_tilt=zeros(length(levs)+1,length(tcbins)-1,length(shearbins)-1);   % The final index in the first dimension will represent the surface center, which remains 0.
  y_tilt=zeros(length(levs)+1,length(tcbins)-1,length(shearbins)-1);   % Indexes 1-6, in this case, will go from 400-->900 hPa.

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
    u=U.u(:,:,:,p);           % U wind
    v=V.v(:,:,:,p);           % V wind
    height=Z.height(:,:,:,p); % Geopotential height

  % .nc SECTION - COMMENT OUT IF USING .mat FORMAT!
%    start=[1 1 1 p]; interval=[length([-10:xgrid:10]) length([-10:ygrid:10]) length(lev) 1];
%    u=ncread(['StormCenteredData/' dset_c '_U.nc'],'u',start,interval);
%    v=ncread(['StormCenteredData/' dset_c '_V.nc'],'v',start,interval);
%    height=ncread(['StormCenteredData/' dset_c '_Z.nc'],'height',start,interval);

  % PREPARE TO ROTATE DATA ACCORDING TO SHEAR VECTOR TO FIND SHEAR-RELATIVE CENTROIDS
    dwnshr=atan2d(shear_v(p),shear_u(p)); lat_box=-10:ygrid:10; lon_box=-10:xgrid:10;
    latbox=lat_box+lat(p); lonbox=lon_box+lon(p);

    if (center_find_method == "VORT") % PREPARE TO CALCULATE VORTICITY FROM U/V WIND
      for z=1:length(levs)       % Do this calculation for each vertical level prescribed near the top.
        vort=zeros(4*x10d+1,4*y10d+1);
        for x=2:4*x10d
          for y=2:4*y10d
            dx=lldistkm([latbox(y) lonbox(x-1)],[latbox(y) lonbox(x+1)]).*1000.0;
            dy=lldistkm([latbox(y-1) lonbox(x)],[latbox(y+1) lonbox(x)]).*1000.0;
            du=u(x,y+1,levs(z))-u(x,y-1,levs(z));
            dv=v(x+1,y,levs(z))-v(x-1,y,levs(z));
            vort(x,y)=(dv./dx)-(du./dy);
          end
        end
        lat_rotate=imrotate(meshgrid(latbox),dwnshr.*(-1)+90,'nearest','crop'); % ROTATE LATITUDE MATRIX
        lon_rotate=imrotate(meshgrid(lonbox),dwnshr.*(-1)+90,'nearest','crop'); % ROTATE LONGITUDE MATRIX
        vo=imrotate(vort,dwnshr.*(-1)+90,'nearest','crop');                     % ROTATE VORTICITY MATRIX

    %   FIND WEIGHTED VORTICITY CENTROID AT EACH LEVEL
        matrix=vo(floor((4*x10d+1)/2)+1-x10d:floor((4*x10d+1)/2)+1+x10d,floor((4*y10d+1)/2)+1-y10d:floor((4*y10d+1)/2)+1+y10d); vortsum=0;
        I=lon_box(floor((4*x10d+1)/2)+1-x10d:floor((4*x10d+1)/2)+1+x10d);
        J=lat_box(floor((4*y10d+1)/2)+1-y10d:floor((4*y10d+1)/2)+1+y10d);
        xsum=0; ysum=0;
        for a=2:length(matrix(:,1))-1
          for b=2:length(matrix(1,:))-1
            if (matrix(a,b) > 0)
              vortsum=vortsum+matrix(a,b);
              xsum=xsum+(I(a).*matrix(a,b));   % Find both the x and y displacement of the vorticity centroid.
              ysum=ysum+(J(b).*matrix(a,b));
            end
          end
        end
        xmax=xsum./vortsum; ymax=ysum./vortsum;
        lon_centroid=lon(p)+xmax;   % Finds exact longitude of centroid
        lat_centroid=lat(p)+ymax;   % Finds exact latitude of centroid

    %   NOW DETERMINE DISTANCES FROM SURFACE CENTER IN X AND Y DIMENSIONS RELATIVE TO SHEAR VECTOR
        [xdist xd2]=lldistkm([lat(p) lon(p)],[lat(p) lon_centroid]);
        if (xmax < 0)
          xdist=xdist.*(-1);
        end
        [ydist yd2]=lldistkm([lat(p) lon(p)],[lat_centroid lon(p)]);
        if (ymax < 0)
          ydist=ydist.*(-1);
        end
        for t=1:length(tcbins)-1
          for s=1:length(shearbins)-1
            if (tcvar(p) >= tcbins(t) && tcvar(p) < tcbins(t+1) && shearvar(p) >= shearbins(s) && shearvar(p) < shearbins(s+1))
              x_tilt(z,t,s)=x_tilt(z,t,s)+xdist;    % We'll add up each snapshot and divide by the total count...
              y_tilt(z,t,s)=y_tilt(z,t,s)+ydist;    % ...to get the composite mean.
              count(z,t,s)=count(z,t,s)+1;
              break;
            end
          end
        end
      end
    elseif (center_find_method == "Z")
      for z=1:length(levs)
    %   NOW FIND THE CENTROID OF MINIMUM GEOPOTENTIAL HEIGHT
        lat_rotate=imrotate(meshgrid(latbox),dwnshr.*(-1)+90,'nearest','crop'); % ROTATE LATITUDE MATRIX
        lon_rotate=imrotate(meshgrid(lonbox),dwnshr.*(-1)+90,'nearest','crop'); % ROTATE LONGITUDE MATRIX
        geopot=height(:,:,levs(z))-mean(mean(height(:,:,levs(z))));
        geopot=geopot.*(-1);   % STRUCTURE SIMILAR TO VORTICITY, ZEROING IN ON "POSITIVE" VALUES.
        vo=imrotate(geopot,dwnshr.*(-1)+90,'nearest','crop');  % ROTATE GEOPOTENTIAL HEIGHT MATRIX

        matrix=vo(floor((4*x10d+1)/2)+1-x10d:floor((4*x10d+1)/2)+1+x10d,floor((4*y10d+1)/2)+1-y10d:floor((4*y10d+1)/2)+1+y10d);
        I=lon_box(floor((4*x10d+1)/2)+1-x10d:floor((4*x10d+1)/2)+1+x10d);
        J=lat_box(floor((4*y10d+1)/2)+1-y10d:floor((4*y10d+1)/2)+1+y10d);
        xsum=0; ysum=0; geosum=0;
        for a=2:length(matrix(:,1))-1
          for b=2:length(matrix(1,:))-1
            if (matrix(a,b) > 0)
              geosum=geosum+matrix(a,b);
              xsum=xsum+(I(a).*matrix(a,b));
              ysum=ysum+(J(b).*matrix(a,b));
            end
          end
        end
        xmax=xsum./geosum; ymax=ysum./geosum;
        lon_centroid=lon(p)+xmax;   % Finds exact longitude of centroid
        lat_centroid=lat(p)+ymax;   % Finds exact latitude of centroid

    %   NOW DETERMINE DISTANCES FROM SURFACE CENTER IN X AND Y DIMENSIONS RELATIVE TO SHEAR VECTOR
        [xdist xd2]=lldistkm([lat(p) lon(p)],[lat(p) lon_centroid]);
        if (xmax < 0)
          xdist=xdist.*(-1);
        end
        [ydist yd2]=lldistkm([lat(p) lon(p)],[lat_centroid lon(p)]);
        if (ymax < 0)
          ydist=ydist.*(-1);
        end
        for t=1:length(tcbins)-1
          for s=1:length(shearbins)-1
            if (tcvar(p) >= tcbins(t) && tcvar(p) < tcbins(t+1) && shearvar(p) >= shearbins(s) && shearvar(p) < shearbins(s+1))
              x_tilt(z,t,s)=x_tilt(z,t,s)+xdist;    % We'll add up each snapshot and divide by the total count...
              y_tilt(z,t,s)=y_tilt(z,t,s)+ydist;    % ...to get the composite mean.
              count(z,t,s)=count(z,t,s)+1;
              break;
            end
          end
        end
      end
    end
    toc
  end

  % NaN out anything that's blank, otherwise compute the composite means for all considered levels above the surface.
  for t=1:length(tcbins)-1
    for s=1:length(shearbins)-1
      for z=1:length(levs)
        if (count(z,t,s) == 0)
          count(z,t,s)=NaN;
          x_tilt(z,t,s)=NaN;
          y_tilt(z,t,s)=NaN;
        else
          x_tilt(z,t,s)=x_bin(z,t,s)./count(z,t,s);
          y_tilt(z,t,s)=y_bin(z,t,s)./count(z,t,s);
        end
      end
    end
  end

  % PLOTTING SECTION - CHANGE DIRECTORY ACCORDINGLY!
  cd CompositePlots
  set(groot,'DefaultFigureColor','white')
  sz=200; sy=['o','^','s','p'];
  cr={[1 0 1],[1 0 0],[0 0 0],[0 0 1]};   % Sets color scheme, marker, and symbol size.

  % We'll do a separate plot for each shear bin, overlaying all TC intensities, intensification rates, etc.
  for s=1:length(shearbins)-1
    h=figure;
    plot(x_tilt(:,1,s),y_tilt(:,1,s),'LineStyle','-','LineWidth',2.0,'Color',cr{1},'Marker',sy(1),'MarkerFaceColor',cr{1},'MarkerEdgeColor',cr{1})
    hold on
    for tcb=2:length(tcbintext)   % Plot as many additional vortex tilt profiles as there are available.
      plot(x_tilt(:,tcb,s),y_tilt(:,tcb,s),'LineStyle','-','LineWidth',2.0,'Color',cr{tcb},'Marker',sy(tcb),'MarkerFaceColor',cr{tcb},'MarkerEdgeColor',cr{tcb})
    end
    hold on
    plot([-105 105],[0 0],'--k','LineWidth',1)  % For now, plot out to 105 km displacement. Can be adjusted depending on exact results.
    plot([0 0],[-105 105],'--k','LineWidth',1)
    hold on
    for c=20:20:100   % Plot circles at 20 km increments to more easily interpret distances.
      viscircles([0 0],c,'Color','k','LineStyle','--','EnhanceVisibility',false,'LineWidth',1);
    end
    xlim([-105 105])
    ylim([-105 105])
    set(gca,'XTick',[-100:20:100])   % Again, may be adjusted as needed.
    set(gca,'YTick',[-100:20:100])
    grid on
    set(gca,'FontSize',14)
    xlabel('Across Shear Distance From Sfc. Center (km)')
    ylabel('Along Shear Distance (km)')
    str1=strcat(dset," Surface-",num2str(lev(levs(1)))," hPa Vortex Tilt");
    str2=shearbintext(s);    % Calls text prescribed above to describe shear magnitude.
    title([str1,str2])
    hold on
    xarrow=0; yarrow=0; uarrow=0; varrow=105;   % Plot shear vector in the positive y direction.
    quiver(xarrow,yarrow,uarrow,varrow,'k','LineWidth',3.5)
    leg=legend(tcbintext,'location','southwest');  % TC bins are listed here. southwest orientation should be out of the way.
    leg.FontSize=11;
    pbaspect([1 1 1])
    hold off
    frameid=sprintf('%02s_VortexTilt_Method_%02s_%02s_%02s_%02d.png',dset,center_find_method,tcplotname,shearplotname,s);
    set(gcf,'inverthardcopy','off')
    print(h,'-dpng',frameid);
    clf
  end
end
