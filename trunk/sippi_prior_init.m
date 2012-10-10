function prior=sippi_prior_init(prior);
% sippi_prior_init Initialize PRIOR structure for SIPPI
%
% Call 
%   prior=sippi_prior_init(prior);
%
% See also sippi_prior
%
disp(sprintf('%s : Initializing prior options',mfilename))


for im=1:length(prior);
    if isfield(prior{im},'Cm');
        if ~isfield(prior{im},'Va');
            prior{im}.Va=prior{im}.Cm;
        end
    end
   
    if ~isfield(prior{im},'type'); 
        disp(sprintf('%s : FATAL ERROR : no ''type'' set for prior %03',mfilename,im));
        return
    end
    prior{im}.type=upper(prior{im}.type);
    
    
    if ~isfield(prior{im},'name'); prior{im}.name=sprintf('%s : m%02d',prior{im}.type,im);end
    
    if ~isfield(prior{im},'x'); prior{im}.x=0;end
    if ~isfield(prior{im},'y'); prior{im}.y=0;end
    if ~isfield(prior{im},'z'); prior{im}.z=0;end
    if ~isfield(prior{im},'dim');
        prior{im}.dim(1)=length(prior{im}.x);
        prior{im}.dim(2)=length(prior{im}.y);
        prior{im}.dim(3)=length(prior{im}.z);
    end
    [prior{im}.xx,prior{im}.yy,prior{im}.zz]=meshgrid(prior{im}.x,prior{im}.y,prior{im}.z);
    
    
    prior{im}.ndim=length(find(prior{im}.dim>1));
    
    %% WRITE TI IF APPLICABLE
    if isfield(prior{im},'TI')
        sgems_write(prior{im}.S.ti_file,prior{im}.TI);
    end
    
    
    
    %% SEQUENTIAL GIBBS OPTIONS
    if ~isfield(prior{im},'seq_gibbs');prior{im}.seq_gibbs.null='';end
    
    %% CHECK FOR GAUSSIAN TYPE PRIOR
    if (strcmp(lower(prior{im}.type),'gaussian'));
        if ~isfield(prior{im}.seq_gibbs,'step_min');prior{im}.seq_gibbs.step_min=0;end
        if ~isfield(prior{im}.seq_gibbs,'step_max');prior{im}.seq_gibbs.step_max=1;end
        if ~isfield(prior{im}.seq_gibbs,'step');prior{im}.seq_gibbs.step=1;end
    end        
    
    %% 
    if ~isfield(prior{im}.seq_gibbs,'type');
        %prior{im}.seq_gibbs.type=1;% BOX RESIM
        prior{im}.seq_gibbs.type=2;% RANDOM POINTS
    end
    if (prior{im}.seq_gibbs.type==1);
        step_min=0;
        try;step_min=prior{im}.x(2)-prior{im}.x(1);end
        try;step_min=min([step_min prior{im}.y(2)-prior{im}.y(1)]);;end
        try;step_min=min([step_min prior{im}.z(2)-prior{im}.z(1)]);;end
        step_max=0;
        try;step_max=max(prior{im}.x)-min(prior{im}.x);;end
        try;step_max=max([step_max max(prior{im}.y)-min(prior{im}.y)]);end
        try;step_max=max([step_max max(prior{im}.z)-min(prior{im}.z)]);end
        if ~isfield(prior{im}.seq_gibbs,'step_min');prior{im}.seq_gibbs.step_min=step_min;end
        if ~isfield(prior{im}.seq_gibbs,'step_max');prior{im}.seq_gibbs.step_max=step_max;end
        if ~isfield(prior{im}.seq_gibbs,'step');prior{im}.seq_gibbs.step=prior{im}.seq_gibbs.step_max;end
    else
        if ~isfield(prior{im}.seq_gibbs,'step_min');prior{im}.seq_gibbs.step_min=1./prod(prior{im}.dim);end
        if ~isfield(prior{im}.seq_gibbs,'step_max');prior{im}.seq_gibbs.step_max=1;end
        if ~isfield(prior{im}.seq_gibbs,'step');prior{im}.seq_gibbs.step=1;end
    end
    
    
    % OPTIONS FOR UPDATING STEP LENTGH
    if ~isfield(prior{im}.seq_gibbs,'i_update_step');
        % Update step length for every i_update_step
        prior{im}.seq_gibbs.i_update_step=50;
    end
    if ~isfield(prior{im}.seq_gibbs,'i_update_step_max');
        % Update step length only for iteration number below i_update_step_max
        prior{im}.seq_gibbs.i_update_step_max=1000;
    end
    
    if ~isfield(prior{im}.seq_gibbs,'n_update_history');
        prior{im}.seq_gibbs.n_update_history=50;
    end
    
    if ~isfield(prior{im}.seq_gibbs,'P_target');
        prior{im}.seq_gibbs.P_target=0.3;
    end
    
    
    %% FFTMA OPTIONS
    if (strcmp(prior{im}.type,'FFTMA'))
        if ~isfield(prior{im},'fft_options');prior{im}.fftma_options.null='';
            if ~isfield(prior{im}.fftma_options,'constant_C');
                prior{im}.fftma_options.constant_C=1;
            end
        end
    end
    
    %% VISIM OPTIONS    
    if (strcmp(prior{im}.type,'VISIM'))
        if ~isfield(prior{im},'V');
            prior{im}.V=visim_init(prior{im}.x,prior{im}.y,prior{im}.z);
        end
        %if isfield(prior{im},'Va');
        %    if isstruct(prior{im}.Va)
        %        prior{im}.V.Va=prior{im}.Va;
        %    else
        %        prior{im}.V.Va=deformat_variogram(prior{im}.Va);
        %    end               
        %end
        if isfield(prior{im},'m0');
           prior{im}.V.m0=prior{im}.m0;
        end        
    end
    
    %% SNESIM OPTIONS    
    if (strcmp(prior{im}.type,'SNESIM'))
        if ~isfield(prior{im},'S');
            prior{im}.S=sgems_get_par('snesim_std');
            prior{im}.S.dim.x=prior{im}.x;            
            prior{im}.S.dim.y=prior{im}.y;
            prior{im}.S.dim.z=prior{im}.z;
        end
        if isfield(prior{im},'ti');
            if isnumeric(prior{im}.ti);
                 prior{im}.S.ti=prior{im}.ti;
            else
                prior{im}.S.ti_file=prior{im}.ti;
            end
            
        else
            prior{im}.ti=prior{im}.S.ti_file;
        end
    end
    
    %% SISIM OPTIONS    
    if (strcmp(prior{im}.type,'SISIM'))
        if ~isfield(prior{im},'S');
            prior{im}.S=sgems_get_par('sisim');
            prior{im}.S.dim.x=prior{im}.x;            
            prior{im}.S.dim.y=prior{im}.y;
            prior{im}.S.dim.z=prior{im}.z;
        end        
    end
    
    %% TARGET DIST
    if isfield(prior{im},'d_target')
        if ~isfield(prior{im},'o_nscore');
            [d_nscore,prior{im}.o_nscore]=nscore(prior{im}.d_target,1,1);
        end
    end
    
    
end

  