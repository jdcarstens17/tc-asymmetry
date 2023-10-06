% FUNCTION SampleSizes - Edited by Jake Carstens 10/6/2023 for manuscript on tropical cyclone asymmetry in reanalyses
%                        Generalized version of script used to generate Figure 1 of manuscript.
% PURPOSE: Use TempestExtremes-derived reanalysis tropical cyclone tracks and intensities, along with wind shear
%          information, to determine composite sample sizes and plot them as bar graphs.
% NOTES: 1. TempestExtremes-derived TC tracks, and wind shear files will be provided
%           in same directory within repository, so filepaths should not need to be altered.
%        2. The user has the choice to use either .mat or .nc format for the wind shear data. TempestExtremes
%           tracks are provided as .txt files.
% PROCEDURE:
% 1. Choose reanalysis dataset (currently ERA5 or CFSR), loading in its TC tracks and associated wind shear.
% 2. Loop through all TC snapshots to place into composite bins based on TC and shear characteristics.
% 3. Plot bar graphs of composite bin sample sizes, with different bar groups referring to the different TC characteristics,
%    and different colors of bars referring to different magnitudes of vertical wind shear.

% FIRST, CHOOSE THE DATASET AND LOAD IN ITS TC TRACKS, INCLUDING STORM ID, YEAR, POSITION, PMIN, AND VMAX.
dset="CFSR"; dset_c='CFSR';
A=readmatrix(['TCTracks/trajectories_' dset_c '.txt']);
id=single(A(:,1)+1); year=single(A(:,2)); lat=A(:,9); pres=A(:,10); wind=A(:,11);
clear A

% LOAD IN WIND SHEAR DATA. WE WILL ONLY NEED THE MAGNITUDE VARIABLE FOR THIS PURPOSE.
% COMMENT OUT ONE OF THE TWO SECTIONS BELOW.
% .mat SECTION
load(['WindShear/Shear_Annulus_200to800km_' dset_c '.mat']) % OTHER OPTIONS FOR SHEAR METHODOLOGY AVAILABLE!
shearvar=shear_magnitude;  % Other variables are shear direction, and zonal/meridional components.
% .nc SECTION
%shearvar=ncread(['WindShear/Shear_Annulus_200to800km_' dset_c '.nc'],'shear_magnitude');

% CHOOSE BINNING SCHEME AND APPROPRIATE PLOTTING LABELS
shearbins=[0:5:10 30];     % Bin values for compositing. Shear is in m/s.
shearbintext=["Low Shear","Mod. Shear","High Shear"];  % This will go in the legend of each bar graph.

% TC CHARACTERISTIC BINNING SCHEMES (WE WILL LOOP THROUGH ALL 4):
% 1 - Intensity - Tropical depression, tropical storm, hurricane, major hurricane
% 2 - Pressure - 10 hPa wide bins starting at < 980 hPa up to > 1000 hPa
% 3 - Intensification Rate - Based on the 6 hour intensity change AFTER the time step in question
% 4 - Latitude - 10 degree bins starting at 0 degrees and going up to 30 (eliminate ET cases)

for tcbinscheme=1:4
  if (tcbinscheme == 1)
    tcbins=[10 17 25 33 100]; tcbintext=["TD","Weak TS","Strong TS","Hurricane"]; % Winds are in m/s.
    tcvar=wind;   % Extracts maximum wind speed from TempestExtremes snapshots.
    tcplotname="VMax";     % This will go in the plot file name.
    xtitle='TC Intensity'; % This will go on the x-axis of the bar graph.
  elseif (tcbinscheme == 2)
    tcbins=[900 980:10:1010]; tcbintext=["< 980 hPa","980-990 hPa","990-1000 hPa","> 1000 hPa"];
    tcvar=pres./100.0;   % Pressure is in Pa in TempestExtremes tracks, so convert to hPa.
    tcplotname="PMin"; xtitle='TC Minimum Pressure';
  elseif (tcbinscheme == 3)
    tcbins=[-10 -2 2 10];    % 2 m/s change over 6 hours equates to about 10 knots in 24 h.
    tcbintext=["Weakening","Steady-State","Intensifying"];
    tcvar=[]; tcplotname="Rate"; xtitle='TC Intensity Change'; % Leave tcvar blank for now since it's set at each step.
  elseif (tcbinscheme == 4)
    tcbins=[0:10:30]; tcbintext=["0-10^o","10-20^o","20-30^o"];
    tcvar=lat; tcplotname="Lat"; xtitle='TC Center Latitude';
  end
  count=zeros(length(tcbins)-1,length(shearbins)-1);  % This simply accumulates snapshots into our composite bins.

% ENTER LOOP TO ASSIGN ALL TC AND SHEAR CHARACTERISTICS INTO BINS, AND PLOT BAR GRAPHS OF THE COMPOSITE SAMPLE SIZE.
  for p=1:length(year)
    if (lat(p) < 0 || lat(p) > 30)   % Neglect extratropical and Southern Hemisphere TCs. Of course, this could be included!
      continue;
    end
    if (tcbinscheme == 3)  % Handle 6-hour intensification rate here, or skip ahead if it's a TC's last timestep
      if (p == length(year))  % Exit the loop if we're at the last snapshot and have no future track.
        break;
      end
      if (id(p+1) ~= id(p))
        continue;
      else
        tcvar(p)=wind(p+1)-wind(p);
      end
    end

%   CHECK TO SEE WHAT COMPOSITE BIN THE SNAPSHOT FALLS INTO, AND ADD THAT SNAPSHOT INTO THE TOTAL COUNT.
    for t=1:length(tcbins)-1
      for s=1:length(shearbins)-1
        if (tcvar(p) >= tcbins(t) && tcvar(p) < tcbins(t+1) && shearvar(p) >= shearbins(s) && shearvar(p) < shearbins(s+1))
          count(t,s)=count(t,s)+1;
          break;   % Break out of this once we find the correct bin.
        end
      end
    end
  end

% PLOTTING SECTION - MAKE SURE TO CHANGE DIRECTORY ACCORDINGLY. NO NEED FOR ADDITIONAL PLOTTING UTILITIES HERE.
  cd Plots
  set(groot,'DefaultFigureColor','white')
  colors=[0 0.6 0     % Low shear: Green
          0.9 0.9 0   % Mod. shear: Yellow
          1 0 0];     % High shear: Red

  h=figure;
  b=bar(count,'FaceColor','flat','EdgeColor','k')  % We'll do this as a bar graph with black outlines, and G/Y/R bars for shear.
  b(1).CData=colors(1,:);
  b(2).CData=colors(2,:);
  b(3).CData=colors(3,:);
  set(gca,'XTick',[1:1:length(tcbins)-1])
  set(gca,'XTickLabel',tcbintext)         % Place the text labels pertaining to each TC bin in the place of numbers on the x-axis.
  set(gca,'FontSize',14)
  xlabel(xtitle)
  ylabel('Total Number of Snapshots')
  str=strcat(dset," Composite Sample Sizes");
  title(str)
  grid on
  leg=legend(shearbintext,'location','best');  % Adjust location if some data are blocked.
  leg.FontSize=11;
  frameid=sprintf('%02s_%02s_SampleSizes.png',dset,tcplotname);
  set(gcf,'inverthardcopy','off')
  print(h,'-dpng',frameid);
  clf
end
