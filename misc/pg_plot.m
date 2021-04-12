% pg_plot: plot plurigaussian transfer function
%
% Call: 
%   [M,x,y]=pg_plot(pg_map,pg_limits);
% See also: pg_transform
%
function [M,x,y]=pg_plot(pg_map,pg_limits,n);

if nargin<2
    pg_limits=[-3 3];
end

if nargin<3
    n=41;
end

[n2,n1]=size(pg_map);

if n2==1;
    x=linspace(pg_limits(1),pg_limits(2),n);
    M=pg_transform(x,pg_map,pg_limits);
    plot(x,M,'.');
    xlabel('Gaussian #1')
    ylabel('index')
    ylim=get(gca,'ylim');
    dy=(ylim(2)-ylim(1)).*.1;
    ylim(1)=ylim(1)-dy;
    ylim(2)=ylim(2)+dy;
    set(gca,'ylim',ylim);
    % 1D
else
    % 2D    
    x=linspace(pg_limits(1),pg_limits(2),n);
    y=linspace(pg_limits(1),pg_limits(2),n);
    [xx,yy]=meshgrid(x,y);
    D{1}=xx;
    D{2}=yy;
    
    DD{1}=xx;DD{2}=yy;
    M=pg_transform(DD,pg_map,pg_limits);
    
    imagesc(x,y,M);
    xlabel('Gaussian #1')
    ylabel('Gaussian #2')
    axis image
    cb=colorbar;
    set(get(cb,'Ylabel'),'String','index')
    grid on
end
%keyboard
