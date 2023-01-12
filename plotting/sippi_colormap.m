% sippi_colormap Default colormap for sippi
%
% Call :
%   sippi_colormap; % the same as sippi_colormap(3);
%
% or :
%   sippi_colormap(1) - Red Green Black
%   sippi_colormap(2) - Red Green Blue Black
%   sippi_colormap(3) - Jet
%   sippi_colormap(4) - Parula
%   sippi_colormap(5) - Geosoft
%
function cmap=sippi_colormap(ic);

if nargin==0
    ic=4;
end
if ic==1
    cmap=(cmap_linear([1 0 0; 0 1 0 ;0 0 0]));
elseif ic==2
    cmap=(cmap_linear([1 0 0; 0 1 0 ; 0 0 1 ;0 0 0]));
elseif ic==3
    cmap=flipud(gray);
elseif ic==4
    if isoctave
      cmap=viridis;
    else 
      cmap=parula;
    end
elseif ic==5
    cmap=cmap_geosoft;
else
    cmap=jet;
end
