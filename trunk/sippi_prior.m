% sippi_prior A priori models for SIPPI 
%
% To generate a realization of the prior model defined by the prior structure use: 
%   [m_propose,prior]=sippi_prior(prior);
%
% To generate a realization of the prior model defined by the prior structure,
% in the vicinity of a current model (using sequential Gibbs sampling) use: 
%   [m_propose,prior]=sippi_prior(prior,m_current);
%
% The following types of a priori models can be used
%   SNESIM  [1D-3D] : based on a multiple point statistical model inferref from a training images. Relies in the SNESIM algorithm
%   SISIM   [1D-3D] : based on Sequential indicator SIMULATION
%   VISIM   [1D-3D] : based on Sequential Gaussian and Direct Sequential simulation
%   FFTMA   [1D-3D] : based on the FFT-MA method (Multivariate Gaussian) 
%   GAUSSIAN   [1D] : 1D generalized gaussian model
%
%
%%%% SIMPLE EXAMPLE %%%
%
%% A simple 2D multivariate Gaissian based prior model based on the 
%% FFT-MA method, can be defined using 
%   id=1;
%   prior{id}.type='FFTMA';
%   prior{id}.name='A SIMPLE PRIOR';
%   prior{id}.x=[0:1:100];
%   prior{id}.y=[0:1:100];
%   prior{id}.m0=10;
%   prior{id}.Va='1 Sph(10)';
%   prior=sippi_prior_init(prior);
%% A realization from this prior model can be generated using
%   m=sippi_prior(prior);
%% This realization can now be plotted using
%   sippi_plot_prior(m,prior);
%% or
%   imagesc(prior{1}.x,prior{1}.y,m{1})
%
%%%% A PRIOR MODEL WITH SEVERAL 'TYPES OF A PRIORI MODEL'
%
%   id=1;
%   prior{id}.type='FFTMA';
%   prior{id}.x=[0:1:100];
%   prior{id}.y=[0:1:100];
%   prior{id}.m0=10;
%   prior{id}.Cm='1 Sph(10)';
%   id=2;
%   prior{id}.type='SISIM';
%   prior{id}.x=[0:1:100];
%   prior{id}.y=[0:1:100];
%   prior{id}.m0=10;
%   prior{id}.Cm='1 Sph(10)';
%   id=3;
%   prior{id}.type='GAUSSIAN';
%   prior{id}.m0=100;
%   prior{id}.std=50;
%   prior{id}.norm=100;
%   prior=sippi_prior_init(prior);
%
%   sippi_plot_model(prior);
%
%%% Sequential Gibbs sampling
%% For more information, see <a href="matlab:web('http://dx.doi.org/10.1007/s10596-011-9271-1')">Hansen, T. M., Cordua, K. S., and Mosegaard, K., 2012. Inverse problems with non-trivial priors - Efficient solution through Sequential Gibbs Sampling. Computational Geosciences</a>.
%
%
% See also: sippi_prior_init, sippi_plot_prior, sippi_prior_set_steplength.m
%
% TMH/2012
%

function [m_propose,prior]=sippi_prior(prior,m_current);

nm=length(prior);
if nargin>1
    m_propose=m_current;
end

%% SELECT WHICH MODEL PARAMETERS TO PERTURB
im_array=[];
for im=1:nm
    if ~isfield(prior{im},'perturb'),prior{im}.perturb=1;end;
    if prior{im}.perturb==1;
        im_array=[im_array im];
    end
end

if nargin==1;
    im_array=1:nm; % SAMPLE ALL MODEL PARAMETERS
end
if isempty(im_array)
    disp(sprintf('%s : no model parameters perturbed...',mfilename))
end

run_fftma=0;
run_snesim=0;

for im=im_array;
    
    if isfield(prior{im},'Cm');
        if ~isfield(prior{im},'Va');
            prior{im}.Va=prior{im}.Cm;
        end
    end
  
    
    if  (strcmp(prior{im}.type,'SISIM'))
        if ~isfield(prior{im},'S');
            prior{im}.S=sgems_get_par('sisim');
        end
        if isfield(prior{im},'marginal_prob');
            prior{im}.S.XML.parameters.Marginal_Probabilities.value=prior{im}.marginal_prob;
            % prior{im}.marginal_prob=ones(1,n_ind)/n_ind;
            n_ind=length(prior{im}.S.XML.parameters.Marginal_Probabilities.value);
            prior{im}.S.XML.parameters.Nb_Indicators.value=n_ind;
        end
        
        
        if isfield(prior{im},'Va');
            va_xml=sgems_variogram_xml(prior{im}.Va);
            prior{im}.S.XML.parameters.Variogram_Median_Ik=va_xml;
        end
        
        % NREAL = 1
        prior{im}.S.XML.parameters.Nb_Realizations.value=1;
        % SEED
        if isfield(prior{im},'seed');
            prior{im}.S.XML.parameters.Seed.value=prior{im}.seed;
        else
            prior{im}.S.XML.parameters.Seed.value=ceil(rand(1).*1e+6);
        end
        
        % REMOVE CONDITIONAL DATA.
        % FIX : NEED TO CHANGE TO HANDLE CONDITIONAL DATA
        if isfield(prior{im}.S,'f_obs')
            prior{im}.S=rmfield(prior{im}.S,'f_obs');
        end
        prior{im}.S.XML.parameters.Hard_Data_Grid.value='';
        prior{im}.S.XML.parameters.Hard_Data_Property.value='';
        
        %% SEQ GIBBS
        if nargin>1
            mgstat_verbose(sprintf('%s : Sequential Gibbs',mfilename),2)
            prior{im}.S=sgems_set_resim_data(prior{im}.S,m_current{im},prior{im}.seq_gibbs.step,prior{im}.seq_gibbs.type);
        end
        
        %% SISIM SIM
        prior{im}.S = sgems_grid(prior{im}.S);
        m_propose{im} = prior{im}.S.D';
        
    elseif  (strcmp(prior{im}.type,'SNESIM')|isfield(prior{im},'S'))
        % SGEMS / SNESIM
        
        % REMOVE CONDITIONAL DATA.
        % FIX : NEED TO CHANGE TO HANDLE CONDITIONAL DATA
        if isfield(prior{im}.S,'f_obs')
            prior{im}.S=rmfield(prior{im}.S,'f_obs');
        end
        prior{im}.S.XML.parameters.Hard_Data.grid='';
        prior{im}.S.XML.parameters.Hard_Data.property='';
        
        %
        prior{im}.S.XML.parameters.Nb_Realizations.value=1;
        
        if isfield(prior{im},'seed');
            prior{im}.S.XML.parameters.Seed.value=prior{im}.seed;
        else
            prior{im}.S.XML.parameters.Seed.value=ceil(rand(1).*1e+6);
        end
        %disp(prior{im}.S.XML.parameters.Seed.value)
        % CHECK FOR SCALING AND ROTATION
        if isfield(prior{im},'rotation')
            prior{im}.S.XML.parameters.Use_Rotation.value=1;
            prior{im}.S.XML.parameters.Use_Global_Rotation.value=1;
            prior{im}.S.XML.parameters.Global_Angle.value=-1*prior{im}.rotation;
        end
        if isfield(prior{im},'scaling')
            if length(prior{im}.scaling)==1; aff=prior{im}.scaling.*[1 1 1]; end
            if length(prior{im}.scaling)==2; aff(3)=1; end
            prior{im}.S.XML.parameters.Use_Affinity.value=1;
            prior{im}.S.XML.parameters.Use_Global_Affinity.value=1;
            prior{im}.S.XML.parameters.Global_Affinity.value=aff;
        end
        
        if nargin>1
            % SEQUENTIAL GIBBS
            if isfield(prior{im},'index_values');
                m = zeros(size(m_current{im}))-1;
                for i=1:length(prior{im}.index_values)
                    try
                        m(find(m_current{im}==prior{im}.m_values(i)))=prior{im}.index_values(i);
                    end
                end
                m_current{im}=m;
            end
            mgstat_verbose(sprintf('%s : Sequential Gibbs',mfilename),2)
            prior{im}.S=sgems_set_resim_data(prior{im}.S,m_current{im},prior{im}.seq_gibbs.step,prior{im}.seq_gibbs.type);
        end
        
        prior{im}.S = sgems_grid(prior{im}.S);
        m_propose{im} = prior{im}.S.D';
        
        if isfield(prior{im},'index_values');
            for i=1:length(prior{im}.index_values)
                m_propose{im}(find(m_propose{im}==prior{im}.index_values(i)))=prior{im}.m_values(i);
            end
        end
        
        
    elseif (strcmp(prior{im}.type,'VISIM')|isfield(prior{im},'V'))
        
        
        % VISIM PRIOR
        
        %% RANDOM SEED
        prior{im}.V.cond_sim=0;
        if isfield(prior{im},'seed');
            prior{im}.V.rseed=prior{im}.seed;
        else
            prior{im}.V.rseed=ceil(rand(1).*1e+6);
        end
        
        if isfield(prior{im},'d_target')
            f_cond=sprintf('d_target_%02d.eas',im);
            if ~exist(f_cond,'file');
                write_eas(f_cond,prior{im}.d_target);
            end
            prior{im}.V.refhist.fname=f_cond;
            prior{im}.V.ccdf=1;
        end
        
        
        %% SEQUENTIAL GIBBS
        if nargin>1
            % SEQUENTIAL GIBBS
            mgstat_verbose(sprintf('%s : Sequential Gibbs',mfilename),2)
            %prior{im}.S=sgems_set_resim_data(prior{im}.S,m_current,prior{im}.seq_gibbs.step,prior{im}.seq_gibbs.type);
            [prior{im}.V, i_resim]=visim_set_resim_data(prior{im}.V,m_current{im},prior{im}.seq_gibbs.step,[],[],prior{im}.seq_gibbs.type);
            if isempty(i_resim)
                prior{im}.V.cond_sim=0;
            else
                prior{im}.V.cond_sim=2;
            end
            
        end
        
        %% RUN VISIM
        prior{im}.V=visim(prior{im}.V);
        m_propose{im} = prior{im}.V.D';
        
        
    elseif (strcmp(prior{im}.type,'FFTMA'))       
        %% FFTMA
        % update VA structure from range/ang/sill if possible
        if ~isstruct(prior{im}.Va); prior{im}.Va=deformat_variogram(prior{im}.Va);end
        
        % UPDATE COVARIANCE PARAMETERS IF THE HAVE BEEN DEFINED
        %[range,rot,sill,Va]=Va2RangeRot(prior{im}.Va);        
        prior{im}.fftma_options.constant_C=1;
        
        % run thorugh the priors that have been updated
        im2_arr=im_array(1:(find(im_array==im)-1));
        for j=im2_arr;
        %for j=1:(im-1);
            
            update_master=0;
            % CHECK IF THIS PRIOR IS A MASTER
            if isfield(prior{j},'prior_master')
                if (prior{j}.prior_master==im);
                    % YES -> current im is master for current j'th prior
                    update_master=1;
                end
            end
            if update_master==1;
                % THIS IM IS A MASTER SO WE CAN UPDATE THE COVARIANCE MODEL
                % IF CHOSEN 
                if strcmp(prior{j}.name,'m0');
                    prior{im}.m0=m_propose{j};
                    prior{im}.fftma_options.constant_C=0;
                elseif strcmp(prior{j}.name,'sill');
                    prior{im}.Va.par1(1)=m_propose{j};
                    prior{im}.fftma_options.constant_C=0;
                elseif strcmp(prior{j}.name,'range_1');
                    range(1)=m_propose{j};
                    prior{im}.Va.par2(1)=range(1);
                    prior{im}.fftma_options.constant_C=0;
                elseif strcmp(prior{j}.name,'range_2');
                    range(2)=m_propose{j};
                    prior{im}.Va.par2(3)=range(2)/range(1);
                    %range(2)=m_propose{j};
                    prior{im}.fftma_options.constant_C=0;
                elseif strcmp(prior{j}.name,'range_3');
                    %range(3)=m_propose{j};
                    prior{im}.fftma_options.constant_C=0;
                elseif strcmp(prior{j}.name,'ang_1');
                    ang(1)=m_propose{j};
                    prior{im}.Va.par2(2)=ang(1);
                    prior{im}.fftma_options.constant_C=0;
                end
            end
        end
        
        % set resim type and step length
        prior{im}.fftma_options.resim_type=prior{im}.seq_gibbs.type;
        prior{im}.fftma_options.lim=prior{im}.seq_gibbs.step;
        
        % IF CONSTANT_C=1, THEN REMOVE C and FFTC
        if prior{im}.fftma_options.constant_C==0;
            try;prior{im}.fftma_options=rmfield(prior{im}.fftma_options,'fftC');end
            try;prior{im}.fftma_options=rmfield(prior{im}.fftma_options,'C');end
        end
        
        
        %if (length(prior{im}.z)==1)&(length(prior{im}.y)==1)
        %    [m_propose{im},prior{im}.fftma_options.z_rand,o]=fft_ma(prior{im}.x,prior{im}.Va,fftma_options);
        %elseif (length(prior{im}.z)==1)
        %    [m_propose{im},prior{im}.fftma_options.z_rand,o]=fft_ma(prior{im}.x,prior{im}.y,prior{im}.Va,fftma_options);
        %else
        %    [m_propose{im},prior{im}.fftma_options.z_rand,o]=fft_ma(prior{im}.x,prior{im}.y,prior{im}.z,fftma_options);
        %end
        [m_propose{im},z_rand,prior{im}.fftma_options]=fft_ma(prior{im}.x,prior{im}.y,prior{im}.z,prior{im}.Va,prior{im}.fftma_options);
        prior{im}.fftma_options.z_rand=z_rand;
        
        
        % PERFORM NORMAL SCORE OF NEEDED        
        if isfield(prior{im},'o_nscore');
            if ~isstruct(prior{im}.Va);
                prior{im}.Va=deformat_variogram(prior{im}.Va);
            end
            Va_par=prior{im}.Va;
            gvar=sum([Va_par.par1]);
            m_propose{im}=m_propose{im}./sqrt(gvar);
            m_propose{im}=inscore(m_propose{im},prior{im}.o_nscore);
        else
            m_propose{im}=m_propose{im}+prior{im}.m0;
        end
        prior{im}.m=m_propose{im};
        
    elseif (strcmp(lower(prior{im}.type),'gaussian'))
        %% 1D GENERALIZED GAUSSIAN
        
        if ~isfield(prior{im},'norm');
            prior{im}.norm=2;
        end
       
        if ~isfield(prior{im},'std');
            if isfield(prior{im},'min')&isfield(prior{im},'max');
                prior{im}.std=(prior{im}.max-prior{im}.min)/2;
            else
                prior{im}.std=1;
            end
        end
        if ~isfield(prior{im},'m0');
            if isfield(prior{im},'min')&isfield(prior{im},'max');
                prior{im}.m0=(prior{im}.max+prior{im}.min)/2;
            else
                prior{im}.m0=0;
            end
        end
        
        if prior{im}.norm==2
            if nargin>1
                gauss_real=(m_current{im}-prior{im}.m0)./prior{im}.std;
                gauss_real_new=randn(1);
                step=pi/2*prior{im}.seq_gibbs.step;
                gauss_real=(cos(step)*gauss_real+sin(step)*gauss_real_new);
                
                %m_propose{im}=gauss_real.*prior{im}.std+prior{im}.m0;
            else
                gauss_real=randn(1);
            end
            m_propose{im}=gauss_real.*prior{im}.std+prior{im}.m0;
        else
            if nargin>1
                if isfield(prior{im},'gauss_real');
                    gauss_real=prior{im}.gauss_real;
                else
                    gauss_real=randn(1);
                end
                gauss_real_new=randn(1);
                %t=sprintf('%5.3f %5.3f',gauss_real,gauss_real_new);
                step=pi/2*prior{im}.seq_gibbs.step;
                gauss_real=(cos(step)*gauss_real+sin(step)*gauss_real_new);
                %t=sprintf('%s %5.3f',t,gauss_real);
                %disp(t);
                
            else
                gauss_real=randn(1);
            end
            
            prior{im}.gauss_real=gauss_real;
            % transform gauss_real to generalized gauss
            g_cdf=normcdf(gauss_real,0,1);
            if ~isfield(prior{im},'ggauss_cdf');
                x=-5:.001:5;nx=length(x);
                pdf=generalized_gaussian(x,0,1,prior{im}.norm,0);
                prior{im}.ggauss_cdf.x=x;
                prior{im}.ggauss_cdf.cdf=cumsum(pdf)./sum(pdf);
            end
            
            % LOCATE VALUE IN GENERALIZED GAUSSIAN DIST.
            jj=find( prior{im}.ggauss_cdf.cdf>0 & prior{im}.ggauss_cdf.cdf<1);
            ggauss_real=interp1(prior{im}.ggauss_cdf.cdf(jj),prior{im}.ggauss_cdf.x(jj),g_cdf);
            if isnan(ggauss_real)
                ii=max(find(prior{im}.ggauss_cdf.cdf<g_cdf));
                if isempty(ii);ii=1;end
                ggauss_real=prior{im}.ggauss_cdf.x(ii);
            end
            
            
            m_propose{im}=prior{im}.m0+prior{im}.std*ggauss_real;
            
        end
        
    else
        disp(sprintf('%s : ''%s'' type prior model not supported',mfilename,prior{im}.type))
    end
    
    
    
    
    %% FIX EXTREME VALUES
    if isfield(prior{im},'min');
        ii=find(m_propose{im}<prior{im}.min);
        m_propose{im}(ii)=prior{im}.min;
        if nargin>2
            m_propose{im}(ii)=m_current{im}(ii);
        end
    end
    if isfield(prior{im},'max');
        ii=find(m_propose{im}>prior{im}.max);
        m_propose{im}(ii)=prior{im}.max;
        if nargin>2
            m_propose{im}(ii)=m_current{im}(ii);
        end
    end
    
    
end


%% UPDATE 'MASTER' PRIOR VARIABLES
%for im=1:nm;
%    if isfield(prior{im},'prior_master');
%        imaster=prior{im}.prior_master;
%        if length(prior)>im
%        prior{imaster}.(prior{im}.name)=m_propose{im};
%        end
%        %disp(sprintf('updated PRIOR%d.%s from PRIOR%d.%s %g',imaster,prior{imaster}.name,im,prior{im}.name,m_propose{im}))
%    end
%end


%% CHECK IF ANY OF THE UPDATES AFFECTS A 'MASTER' PRIOR VARIABLE
for im=im_array;
    if isfield(prior{im},'prior_master');
        imaster=prior{im}.prior_master;
        if (length(prior)<=imaster)
            prior{imaster}.run=1;
        end
    end
end

%% REMOVE INFORMATION OF MASTER RUN IF SET
%for i=1:nm
%    if (isfield(prior{im},'run'))
%        prior{imaster}=rmfield(prior{imaster},'run')
%    end
%end

% %%
% %% 'UPDATE MASTER' 
% %% FOR NOW ONLY FFTMA CAN BE USED AS A MASTER; 
% %%
for im=1:nm;

    
    if (isfield(prior{im},'run'))
        %disp(sprintf('RUN MASTER %d (%s) AGAIN',im,prior{im}.name));
        if strcmp(prior{im}.type,'FFTMA')
            
            if ~isstruct(prior{im}.Va); prior{im}.Va=deformat_variogram(prior{im}.Va);end
            
            % update covariance from child priors if any
            for j=1:(im-1);
                update_master=0;
                % CHECK IF THIS PRIOR IS A MASTER
                if isfield(prior{j},'prior_master')
                    if (prior{j}.prior_master==im);
                        % YES -> current im is master for current j'th prior
                        update_master=1;
                    end
                end
                if update_master==1;
                    % THIS IM IS A MASTER SO WE CAN UPDATE THE COVARIANCE MODEL
                    % IF CHOSEN
                    if strcmp(prior{j}.name,'m0');
                        prior{im}.m0=m_propose{j};
                        prior{im}.fftma_options.constant_C=0;
                    elseif strcmp(prior{j}.name,'sill');
                        prior{im}.Va.par1(1)=m_propose{j};
                        prior{im}.fftma_options.constant_C=0;
                    elseif strcmp(prior{j}.name,'range_1');
                        range(1)=m_propose{j};
                        prior{im}.Va.par2(1)=range(1);
                        prior{im}.fftma_options.constant_C=0;
                    elseif strcmp(prior{j}.name,'range_2');
                        range(2)=m_propose{j};
                        prior{im}.Va.par2(3)=range(2)/range(1);
                        %range(2)=m_propose{j};
                        prior{im}.fftma_options.constant_C=0;
                    elseif strcmp(prior{j}.name,'range_3');
                        %range(3)=m_propose{j};
                        prior{im}.fftma_options.constant_C=0;
                    elseif strcmp(prior{j}.name,'ang_1');
                        ang(1)=m_propose{j};
                        prior{im}.Va.par2(2)=ang(1);
                        prior{im}.fftma_options.constant_C=0;
                    end
                end
            end
            
            % set resim type and step length
            prior{im}.fftma_options.resim_type=prior{im}.seq_gibbs.type;
            % RANDOM NUMBERS HAVE ALLREADY BEEN UPDATED IF THE FFTMA MASTER
            % PRIOR HAS BEEN PERTURBEN. HERE ONLY THE COVARIANCE 
            % PROPERTIES ARE PERTURBED! (random numbers are fixed)
            prior{im}.fftma_options.lim=0; 
            % IF CONSTANT_C=1, THEN REMOVE C and FFTC
            if prior{im}.fftma_options.constant_C==0;
                try;prior{im}.fftma_options=rmfield(prior{im}.fftma_options,'fftC');end
                try;prior{im}.fftma_options=rmfield(prior{im}.fftma_options,'C');end
            end
            
            
            %if (length(prior{im}.z)==1)&(length(prior{im}.y)==1)
            %    [m_propose{im},prior{im}.fftma_options.z_rand,o]=fft_ma(prior{im}.x,prior{im}.Va,fftma_options);
            %elseif (length(prior{im}.z)==1)
            %    [m_propose{im},prior{im}.fftma_options.z_rand,o]=fft_ma(prior{im}.x,prior{im}.y,prior{im}.Va,fftma_options);
            %else
            %    [m_propose{im},prior{im}.fftma_options.z_rand,o]=fft_ma(prior{im}.x,prior{im}.y,prior{im}.z,fftma_options);
            %end
            [m_propose{im},z_rand,prior{im}.fftma_options]=fft_ma(prior{im}.x,prior{im}.y,prior{im}.z,prior{im}.Va,prior{im}.fftma_options);
            prior{im}.fftma_options.z_rand=z_rand;
        
        
            % PERFORM NORMAL SCORE OF NEEDED
            if isfield(prior{im},'o_nscore');
                if ~isstruct(prior{im}.Va);
                    prior{im}.Va=deformat_variogram(prior{im}.Va);
                end
                Va_par=prior{im}.Va;
                gvar=sum([Va_par.par1]);
                m_propose{im}=m_propose{im}./sqrt(gvar);
                m_propose{im}=inscore(m_propose{im},prior{im}.o_nscore);
            else
                m_propose{im}=m_propose{im}+prior{im}.m0;
            end
            prior{im}.m=m_propose{im};
          
        end
    end
end

%% REMOVE INFORMATION OF MASTER RUN IF SET
for im=1:nm
    if (isfield(prior{im},'run'))
        prior{im}=rmfield(prior{imaster},'run');
    end
end

