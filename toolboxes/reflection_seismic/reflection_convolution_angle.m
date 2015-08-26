% reflection_convolution_angle : convolution of vp,vs,rho profiles
% using full zoeppritz, Shuey
function [seis]=reflection_convolution_angle(vp,vs,rho,angle,wavelet,type)
if nargin<6
    type='full_zoeppritz';
    type='shuey';
    type='shuey_2_term';
end
ntraces=size(vp,2);
ndata=size(vp,1);
nangle=length(angle);

seis=zeros(ndata*nangle,ntraces);

if strcmp(lower(type),'full_zoeppritz');
    %% FULL ZOEPPRITZ
    for i=1:ntraces;
        % OBTAIN RPP SERIES
        rpp=zeros(ndata,nangle);
        for j=2:ndata
            rt_ps=zoepprits_forward(vp(j-1,i),vs(j-1,i),rho(j-1,i),vp(j,i),vs(j,i),rho(j,i),angle)';
            rpp(j,:)=rt_ps(1,:);
        end
        
        % CONVOLUTION
        seis_data=reflectivity_convolution(rpp,wavelet);
        seis(:,i)=seis_data(:);
        
    end
elseif strfind(lower(type),'shuey');
    %% LINEAR AKI AND RICHARDS
    
    
    d_vp=diff(vp);
    d_vs=diff(vs);
    d_rho=diff(rho);
    vp_av=(vp(2:end,:)+vp(1:end-1,:))/2;
    vs_av=(vs(2:end,:)+vs(1:end-1,:))/2;
    rho_av=(rho(2:end,:)+rho(1:end-1,:))/2;
    
    R=.5*(d_vp./vp_av + d_rho./rho_av);
    G=.5*(d_vp./vp_av) -2*(vs_av.^2./vp_av.^2).*(  d_rho./rho_av + 2*d_vs./vs_av );
    F=.5*(d_vp./vp_av);
    X=2*vs_av.^2./vp_av;
    arad=pi*angle./180;
    
    for i=1:ntraces
        rpp=zeros(ndata,nangle);
        for k=1:nangle
            if strfind(lower(type),'shuey_2_term');
                rpp(2:end,k)=R(:,i) + G(:,i)*sin(arad(k)).^2;
            else
                rpp(2:end,k)=R(:,i) + G(:,i)*sin(arad(k)).^2 + F(:,i)*(tan(arad(k)).^2 - sin(arad(k)).^2);
            end
        end
        seis_data=reflectivity_convolution(rpp,wavelet);
        seis(:,i)=seis_data(:);
    end
    
elseif strcmp(lower(type),'bulandomre');
    
    
    
end

%     for j=1:ntraces;
%          % OBTAIN RPP SERIES
%         rpp_ang=zeros(ndata,nangle);
%         for j=2:ndata
%             rt_ps=zoepprits_forward(vp(j-1,i),vs(j-1,i),rho(j-1,i),vp(j,i),vs(j,i),rho(j,i),angle)';
%             rpp(j,:)=rt_ps(1,:);
%         end
%
%         % CONVOLUTION
%         seis_data=reflection_coefficient_to_anglegather(rpp,wavelet);
%         seis(:,i)=seis_data(:);
%
%     end
%


end
