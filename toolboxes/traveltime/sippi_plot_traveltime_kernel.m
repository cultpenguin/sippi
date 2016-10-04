% sippi_plot_traveltime_kernel: plot the forward kernel (if it exists) on
% top of a realization of the prior
%
% Call:
%    sippi_plot_traveltime_kernel(forward,prior);
%    sippi_plot_traveltime_kernel(forward,prior,m);

function sippi_plot_traveltime_kernel(forward,prior,m);

if nargin<2
    disp(sprintf('%s: needs at least 2 inputs',mfilename));
    help(mfilename);
    return
end

fn=get(gcf,'Number');
if nargin<3;
    m=sippi_prior(prior);
end
sippi_plot_prior(prior,m,1,1,fn);
str=get(get(gca,'title'),'string');
title([str,' - Source-Receiver locations'])
hold on
if isfield(forward,'sources');
    plot(forward.sources(:,1),forward.sources(:,2),'r*')
    plot(forward.receivers(:,1),forward.receivers(:,2),'ko')
    plot([forward.sources(:,1) forward.receivers(:,1)]',[forward.sources(:,2) forward.receivers(:,2)]','k-')
end
hold off

if isfield(forward,'G');
    figure_focus(fn+1);
    subplot(1,2,1);
    imagesc(prior{1}.x,prior{1}.y,reshape(forward.G(1,:), length(prior{1}.y), length(prior{1}.x) ))
    xlabel('X');ylabel('Y')
    axis image;
    title('kernel 1')
    subplot(1,2,2);
    imagesc(prior{1}.x,prior{1}.y,reshape(sum(forward.G), length(prior{1}.y), length(prior{1}.x) ))
    axis image
    xlabel('X');ylabel('Y')
    title('sum of all kernels')
end


