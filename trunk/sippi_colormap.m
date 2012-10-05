% sippi_colormap Default colormap for sippi
%
% Call :
%   sippi_colormap; % the same as sippi_colormap(3);
%
% or :
%   sippi_colormap(1) - Red Green Black
%   sippi_colormap(2) - Red Green Blue Black
%   sippi_colormap(3) - Jet
%
function cmap=sippi_colormap(ic);

if nargin==0
    ic=3;
end
if ic==1
    cmap=(cmap_linear([1 0 0; 0 1 0 ;0 0 0]));
elseif ic==2
    cmap=(cmap_linear([1 0 0; 0 1 0 ; 0 0 1 ;0 0 0]));
else
    cmap=jet;
end