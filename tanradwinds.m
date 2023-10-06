function [utan,urad]=tanradwinds(u,v)
% tanradwinds.m - Written by Jake Carstens 9/7/2022
% Takes u and v wind components and converts to tangential and radial winds.
% Data must be in 2 dimensions, i.e. you would need to make separate calls to this
% function to calculate a vertical profile of radial/tangential winds.

  [nx,ny]=size(u); xc=floor(nx/2)+1; yc=floor(ny/2)+1; % Establish grid size/TC center
  [xx,yy] = meshgrid(1:nx,1:ny); theta = atan2(xx-xc,yy-yc);

% Calculating radial and tangential winds on cartesian grid
  urad = u.*cos(theta)+v.*sin(theta); % radial wind
  utan = -u.*sin(theta)+v.*cos(theta); % tangential wind
end
