% unknown_pleasures: Simulaing the cover photo from Joy Divisions Unknown
% Pleasures album, suing SIPPI.
% cover from unknown_pleasures
clear all

%% Setup two prior structures
x=linspace(0.1,.9,101);
d_target=rand(1,1000);

nl=60; % number of lines

ip=1;
prior{ip}.type='FFTMA';
prior{ip}.x=x;;
prior{ip}.y=1:1:nl;
prior{ip}.m0=0;
prior{ip}.Va='0.1 Sph(.1,90,0.001)';
prior{ip}.d_target=d_target/20;

ip=2;
prior{ip}.type='FFTMA';
prior{ip}.x=x;;
prior{ip}.y=1:1:nl;
prior{ip}.m0=0;
prior{ip}.Va='.99 Gau(.06,90,0.001) + .01 Sph(.1,90,0.001)';
prior{ip}.d_target=3*d_target;


%% Compute a padding matrix, such that prior{2} is only used 'in the middle' 
pad=zeros(size(x));


x1=0.3;ix1=find(x>=x1);ix1=min(ix1);
x2=0.4;;ix2=find(x>=x2);ix2=ix2(1);
xx1=ix1:1:ix2;
fadein=sin(interp1([ix1 ix2],[0 pi/2],xx1));
pad(find(x>x1&x<(1-x2)))=1;
pad(xx1)=fadein;

xx2=fliplr(length(x)-xx1);
pad(xx2)=fliplr(fadein);

% add a small fading from left to right
linpad=linspace(1.3,0.0,length(x));
pad=pad.*linpad;

pad=repmat(pad,nl,1);

%% Generate first realization.
[m,prior]=sippi_prior(prior);
mm=m{1}+m{2}.*pad;


%% simulate 'unknown pleasures' cover

ax=[-1 2 -20 nl+20];
prior{1}.seq_gibbs.step=0.1;
prior{2}.seq_gibbs.step=0.01;
prior{1}.fftma_options.constant_C=0;
prior{2}.fftma_options.constant_C=0;

save_movie=0;
save_png=0;

if (save_movie==1);
    vidObj = VideoWriter('unknown_pleasures','MPEG-4');
    vidObj.FrameRate = 30;
    vidObj.Quality = 95;
    open(vidObj)
end

nsim=800;
figure(1);
subfigure(2,2,1)

for isim=1:nsim;

    % setting aniso
    aniso=max([40 400-(isim/(nsim/4))*400]);
    prior{1}.Va(end).par2(3)=aniso;
    prior{2}.Va(1).par2(3)=aniso;
    % setting variance
    scale=min([1 isim/(nsim/4)]);
    %prior{2}.Va(1).par1=scale*.99;
    
    disp(sprintf('isim=%03d/%03d, aniso=%g, scale=%g',isim,nsim,aniso,scale));

    [m,prior]=sippi_prior(prior,m);
    mm=m{1}+scale*m{2}.*pad;
    
    figure_focus(1);clf
    colordef black;
    fill([ax(1) ax(2) ax(2) ax(1)],[ax(3) ax(3) ax(4) ax(4)],[0 0 0]);
    hold on
    
    
    for i=nl:-1:1;
        x_fill=[x,max(x),min(x),min(x)];
        d=i+mm(i,:);
        d_fill=[d,0,0,d(1)];
        f{i}=fill(x_fill,d_fill,'black');
        set(f{i},'linewidth',.1)
        hold on
        p=plot(x_fill,d_fill,'r-');        
    end
    plot([0.1 0.1 0.9 0.9],[nl 0 0 nl],'k-')
    hold off
    axis(ax);
    axis off
    drawnow;
    if (save_png==1)
        print_mul(sprintf('unknown_pleasures_%3d',isim),1)
    end
    if (save_movie==1);
        currFrame = getframe(gcf);
        writeVideo(vidObj,currFrame);
    end
end
if (save_movie==1);
    try
        close(vidObj)
    end
end
