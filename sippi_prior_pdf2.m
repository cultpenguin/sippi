% sippi_prior_pdf: sampled from 2D PDF
function [m_propose,prior]=sippi_prior_pdf2(prior,m_current,ip);


if nargin == 0;
    mu = [6 12;-3 5];
    sigma = cat(3,[4 15],[5 1]); % 1-by-2-by-2 array
    mu = [-5 5;5 16];
    sigma = cat(3,[3 4],[3 14]); % 1-by-2-by-2 array
    gm = gmdistribution(mu,sigma);
    
    nx=131;
    ny=141;
    x=linspace(-10,11,nx);
    y=linspace(0,20,ny);
    [xx,yy]=meshgrid(x,y);
    P=reshape(pdf(gm,[xx(:),yy(:)]),ny,nx);
    prior{1}.type='pdf2';
    prior{1}.x=[1,2]; % set the correct 'dimension' / fix automatic
    prior{1}.pdf=P;
    prior{1}.pdf_x=x;
    prior{1}.pdf_y=y;
    
    [m,prior]=sippi_prior(prior);

    prior{1}.seq_gibbs.step=0.05;
    N=1000;
    disp(sprintf('%s: Generating %d realizations.',mfilename,N))
    sim=ones(N,2).*NaN;
    for i=1:N;
        [m,prior]=sippi_prior(prior,m);
        sim(i,:)=m{1};
    end

    for i=1:1:N;
        imagesc(prior{1}.pdf_x,prior{1}.pdf_y,P);
        hold on;
        plot(sim(1:i,1),sim(1:i,2),'w.','MarkerSize',6);
        hold off;
        drawnow;
    end
    
    m_propose = m;
    return
end

if nargin<3;
    ip=1;
end

if ~isfield(prior{ip},'init')
    prior=sippi_prior_init(prior);
end

if length(prior{ip}.x)==1
    prior{ip}.x=[1,2];
end


if ~isfield(prior{ip},'pdf')
end

if ~isfield(prior{ip},'pdf_x')
    nx=size(prior{ip}.pdf,2);
    prior{ip}.pdf_x=1:1:nx;
end

if ~isfield(prior{ip},'pdf_y')
    ny=size(prior{ip}.pdf,1);
    prior{ip}.pdf_y=1:1:ny;
end

if ~isfield(prior{ip},'p')
    prior{ip}.p=[rand(1),rand(1)];
end

perturb=1;
if nargin==1;
    perturb=0;
else
    if prior{ip}.seq_gibbs.step==1
        perturb=0;
    end
end


if perturb==0
    prior{ip}.p=[rand(1),rand(1)];
else
    p_org=prior{ip}.p;
    p=p_org+randn(1,2).*prior{ip}.seq_gibbs.step;
    p(find(p<0))=p_org(find(p<0));
    p(find(p>1))=p_org(find(p>1));
    prior{ip}.p=p;
end


% XDIM
CPDF_x=cumsum(sum(prior{ip}.pdf,1));CPDF_x=CPDF_x./max(CPDF_x(:));
ix_arr=find(CPDF_x>=prior{ip}.p(1));
ix=ix_arr(1);
ix1=max([1 ix-1]);
ix2=min([length(prior{1}.pdf_x) ix+1]);
ixx=ix1:ix2;
try
    x_sim=interp1(CPDF_x(ixx),prior{1}.pdf_x(ixx),prior{1}.p(1));
    if isnan(x_sim)
        x_sim=interp1(CPDF_x(ixx),prior{1}.pdf_x(ixx),prior{1}.p(1),'nearest','extrap');
    end    

    %x_sim=interp1(CPDF_x,x,p(1));
catch
    x_sim=prior{1}.pdf_x(ix);
end

% YDIM
CPDF_y=cumsum(sum(prior{1}.pdf(:,ix),2));CPDF_y=CPDF_y./max(CPDF_y(:));
iy_arr=find(CPDF_y>=prior{1}.p(2));
iy=iy_arr(1);
iy1=max([1 iy-1]);
iy2=min([length(prior{1}.pdf_y) iy+1]);
iyy=iy1:iy2;
try
    y_sim=interp1(CPDF_y(iyy),prior{1}.pdf_y(iyy),prior{1}.p(2));
    if isnan(y_sim)
        y_sim=interp1(CPDF_y(iyy),prior{1}.pdf_y(iyy),prior{1}.p(2),'nearest','extrap');
    end    
catch
    y_sim=prior{1}.pdf_y(iy);
end

m_propose{1}=[x_sim;y_sim];

