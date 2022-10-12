% sippi_plot_traveltime_kernel: plot the forward kernel (if it exists) on
% top of a realization of the prior
%
% Call:
%    sippi_plot_traveltime_kernel(forward,prior);
%    sippi_plot_traveltime_kernel(forward,prior,m);

function sippi_plot_traveltime_kernel(forward,prior,m,pl_kernel,i_use);

if nargin<2
    disp(sprintf('%s: needs at least 2 inputs',mfilename));
    help(mfilename);
    return
end

fn=get(gcf,'Number');
if nargin<3;
    m=sippi_prior(prior);
end
if nargin<4;
    pl_kernel=1;
end
if nargin<5;
    i_use = 1:size(forward.sources,1);
end


sippi_plot_prior(prior,m,1,1,fn);
str=get(get(gca,'title'),'string');
title([str,' - Source-Receiver locations'])
hold on
if isfield(forward,'sources');
    plot(forward.receivers(i_use,1),forward.receivers(i_use,2),'k.','MarkerSize',18)
    plot(forward.sources(i_use,1),forward.sources(i_use,2),'r*')
    plot([forward.sources(i_use,1) forward.receivers(i_use,1)]',[forward.sources(i_use,2) forward.receivers(i_use,2)]','k-','LineWidth',.1)
end
hold off


if isfield(forward,'G')&&(pl_kernel==1);
    figure_focus(fn+1);
    subplot(1,2,1);
    imagesc(prior{1}.x,prior{1}.y,reshape(forward.G(i_use(1),:), length(prior{1}.y), length(prior{1}.x) ))
    hold on
    plot(forward.sources(1,1),forward.sources(1,2),'r*')
    plot(forward.receivers(1,1),forward.receivers(1,2),'ko')
    hold off
    xlabel('X');ylabel('Y')
    axis image;
    title('kernel 1')
    subplot(1,2,2);
    imagesc(prior{1}.x,prior{1}.y,reshape(sum(forward.G(i_use,:)), length(prior{1}.y), length(prior{1}.x) ))
    axis image
    xlabel('X');ylabel('Y')
    title('sum of all kernels')
end


