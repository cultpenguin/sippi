% pg_plot: plot plurigaussian transfer function
%
% See also: pg_transform
%
function [M,x1,y1]=pg_plot(pg_map,pg_limits);

if nargin<2
    pg_limits=[-3 3];
end

[n2,n1]=size(pg_map);

x1=linspace(pg_limits(1),pg_limits(2),n1);
y1=linspace(pg_limits(1),pg_limits(2),n2);
[xx1,yy1]=meshgrid(x1,y1);

lim=linspace(pg_limits(1),pg_limits(2),1001);
[xx,yy]=meshgrid(lim,lim);

if n2==1;
    M=interp1(x1,pg_map,lim,'next');
    plot(lim,M,'.');
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
    M=griddata(xx1(:),yy1(:),pg_map(:),xx,yy,'nearest');

    imagesc(x1,y1,M);
    xlabel('Gaussian #1')
    ylabel('Gaussian #2')
    axis image
    cb=colorbar;
    set(get(cb,'Ylabel'),'String','index')
    grid on
end
%keyboard
