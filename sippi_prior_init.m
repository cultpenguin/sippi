function prior=sippi_prior_init(prior,im_array);
% sippi_prior_init Initialize PRIOR structure for SIPPI
%
% Call
%   prior=sippi_prior_init(prior);
%
% See also sippi_prior
%

if nargin<2
    im_array=1:length(prior);
end

for im=im_array
    % only perform initialization of not allready done.
    if isfield(prior{im},'init')
        if prior{im}.init==1
            im_array=setxor(im_array,im);
        end
    end
end

for im=im_array
    
    
    prior{im}.init=0;
    
    %% AUTOMATICALLY CHECK FOR CM VERSIS VA FIELD
    %% TMH: CM should be deafult..
    if isfield(prior{im},'Cm');
        if ~isfield(prior{im},'Va');
            prior{im}.Va=prior{im}.Cm;
        end
    end
    if isfield(prior{im},'Va');
        if ~isfield(prior{im},'m0');
            prior{im}.m0=0;
        end
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CHECK FOR TYPE    
    if ~isfield(prior{im},'type');
        sippi_verbose(sprintf('%s : FATAL ERROR : no ''type'' set for prior %03',mfilename,im));
        break
    end
    prior{im}.type=upper(prior{im}.type);
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SET OPTIONAL NAME FIELD IF NOT SET
    if ~isfield(prior{im},'name'); prior{im}.name=sprintf('%s : m%02d',prior{im}.type,im);end
    sippi_verbose(sprintf('%s : setting name for prior{%d} as ''%s''',mfilename,im,prior{im}.name),2);
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % SET DIMENSIONS OF PRIOR IF NOT ALLREADY SET
    if ~isfield(prior{im},'x'); prior{im}.x=0;end
    if ~isfield(prior{im},'y'); prior{im}.y=0;end
    if ~isfield(prior{im},'z'); prior{im}.z=0;end
    if ~isfield(prior{im},'dim');
        prior{im}.dim(1)=length(prior{im}.x);
        prior{im}.dim(2)=length(prior{im}.y);
        prior{im}.dim(3)=length(prior{im}.z);
    end
    [prior{im}.xx,prior{im}.yy,prior{im}.zz]=meshgrid(prior{im}.x,prior{im}.y,prior{im}.z);

    % NDATA FOR EACH DIM
    try,prior{im}.lim(1)=max(prior{im}.x)-min(prior{im}.x);end
    try,prior{im}.lim(2)=max(prior{im}.y)-min(prior{im}.y);end
    try,prior{im}.lim(3)=max(prior{im}.z)-min(prior{im}.z);end
    
    prior{im}.ndim=length(find(prior{im}.dim>1));
    
    if (~isfield(prior{im},'cax'));
        if ((isfield(prior{im},'min'))&&(isfield(prior{im},'max')));
            prior{im}.cax=[prior{im}.min prior{im}.max];
        end
        
    end
    
    
    
    %%
    % SEQUENTIAL GIBBS OPTIONS
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~isfield(prior{im},'seq_gibbs');prior{im}.seq_gibbs.null='';end
   
    if (strcmp(lower(prior{im}.type),'gaussian'));
        if ~isfield(prior{im}.seq_gibbs,'step_min');prior{im}.seq_gibbs.step_min=0;end
        if ~isfield(prior{im}.seq_gibbs,'step_max');prior{im}.seq_gibbs.step_max=1;end
        if ~isfield(prior{im}.seq_gibbs,'step');prior{im}.seq_gibbs.step=1;end
    end

    
    %
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
        if ~isfield(prior{im}.seq_gibbs,'step_min');prior{im}.seq_gibbs.step_min=1./(prod(prior{im}.dim*2));end
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
    if (strcmp(upper(prior{im}.type),'FFTMA'))
        if ~isfield(prior{im},'fft_options');prior{im}.fftma_options.null='';
            if ~isfield(prior{im}.fftma_options,'constant_C');
                prior{im}.fftma_options.constant_C=1;
            end
        end
    end

    %%
    % BELOW HERE COMES VARIOUS SANITY CHECKS DIFFERENT PRIOR TYPES
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %% UNIFORM TYPE PRIOR SANITY CHECK
    if (strcmp(upper(prior{im}.type),'UNIFORM'))
        if ~isfield(prior{im},'min'); 
            prior{im}.min=0;
            sippi_verbose(sprintf('%s : setting prior{%d}.min=%g (uniform minimal value)',mfilename,im,prior{im}.min ));
        end
        if ~isfield(prior{im},'max'); 
            prior{im}.max=prior{im}.min+1;
            sippi_verbose(sprintf('%s : setting prior{%d}.max=%g (uniform maximum value)',mfilename,im,prior{im}.max ));
        end
    end
    
    %% GAUSSIAN TYPE PRIOR SANITY CHECK
    if (strcmp(lower(prior{im}.type),'gaussian'));
        % m0
        if ~isfield(prior{im},'m0')&&(~isfield(prior{im},'d_target'))
            if isfield(prior{im},'min')&isfield(prior{im},'max');
                prior{im}.m0=(prior{im}.max+prior{im}.min)/2;
            else
                prior{im}.m0=0;
            end
            sippi_verbose(sprintf('%s : setting prior{%d}.m0=%g (prior mean)',mfilename,im,prior{im}.m0 ));
        end
        
        % std
        if ~isfield(prior{im},'std')&&(~isfield(prior{im},'d_target'))
            if isfield(prior{im},'min')&isfield(prior{im},'max');
                prior{im}.std=(prior{im}.max-prior{im}.min)/2.300;
            else
                prior{im}.std=1;
            end
            sippi_verbose(sprintf('%s : setting prior{%d}.std=%g (prior standard deviation)',mfilename,im,prior{im}.std ));
        end
    end
    
    
    %% VISIM OPTIONS
    if (strcmp(upper(prior{im}.type),'VISIM'))
        visim_clean;
       
        if ~isfield(prior{im},'visim_id');
            prior{im}.visim_id=im;
        end    
       
        if ~isfield(prior{im},'V');
            prior{im}.V=visim_init(prior{im}.x,prior{im}.y,prior{im}.z);
        end
        if ~isfield(prior{im},'m0');
            prior{im}.m0=0;
            sippi_verbose(sprintf('%s : setting prior{%d}.m0=%g (prior mean)',mfilename,im,prior{im}.m0 ));
        end
        if isfield(prior{im},'m0');
            if length(prior{im}.m0)>1
                prior{im}.V.gmean=prior{im}.m0(1);
                txt=sprintf('%s : VISIM only supprt a constant mean\n',mfilename);
                txt=sprintf('%s :   setting prior{%d}.m0=',im,prior{im}.m0);
                sippi_verbose(txt);
            end
            prior{im}.V.gmean=prior{im}.m0;
        end
        try 
             f_cond=sprintf('d_target_%02d.eas',im);
             if exist([pwd,filesep,f_cond],'file')
                 delete(f_cond);
             end
        end
    end
    
    %% CHOLESKY
    if (strcmp(upper(prior{im}.type),'CHOLESKY'))
        if ~isfield(prior{im},'Cm')&~isfield(prior{im},'Va')&~isfield(prior{im},'Cmat')
            prior{im}.Cm='1 Sph(1)';
            txt=sprintf('%s : No covariance model set, using prior{%d}.Cm=''%s''',mfilename,im,prior{im}.Cm);
            sippi_verbose(txt);
        end
        
        % check for correct size of Cmat
        if isfield(prior{im},'Cmat');
            if (prod(prior{im}.dim)~=size(prior{im}.Cmat,1))|(prod(prior{im}.dim)~=size(prior{im}.Cmat,2))
                txt=sprintf('%s : FATAL ERROR: prior{%d}.Cmat is of size [%d,%d],\n',mfilename,im,size(prior{im}.Cmat,2),size(prior{im}.Cmat,1));
                txt=sprintf('%s%s   but should be of size  [%d,%d]',txt,char(32*ones(1,length(mfilename))),prod(prior{im}.dim),prod(prior{im}.dim));
                sippi_verbose(txt);
            end
        end
        
    end
    

    
    %% FFTMA
    if (strcmp(upper(prior{im}.type),'FFTMA'))
        if ~isfield(prior{im},'Cm')&~isfield(prior{im},'Va');
            prior{im}.Cm='1 Sph(1)';
            txt=sprintf('%s : No covariance model set, using prior{%d}.Cm=''%s''',mfilename,im,prior{im}.Cm);
            sippi_verbose(txt);
        end       
    end
        
    %% FFTMA and CHOLESKY
    if (strcmp(upper(prior{im}.type),'FFTMA'))||(strcmp(upper(prior{im}.type),'CHOLESKY'))
        % THE CHECKS BELOW ARE SIMILAT FOR FFTMA AND CHOLESKY TYPE PRIORS
        % PERHAPS CHECK BOTH AT THE SAME TIME
        if ~isfield(prior{im},'m0')
            prior{im}.m0=0;
            txt=sprintf('%s : No mean set, using prior{%d}.m0=%g',mfilename,im,prior{im}.m0);
            sippi_verbose(txt);
        end
        % Check for correct size of m0
        if isfield(prior{im},'m0');
            if prod(size(prior{im}.m0))~=1
                if (size(prior{im}.m0,2)~=prior{im}.dim(1))|(size(prior{im}.m0,1)~=prior{im}.dim(2))|(size(prior{im}.m0,3)~=prior{im}.dim(3))
                    txt=sprintf('%s : FATAL ERROR: prior{%d}.m0, should be either a scalar\n',mfilename,im);
                    txt=sprintf('%s%s   or of size=[%g,%g,%g]',txt,char(32*ones(1,length(mfilename))),prior{im}.dim(2),prior{im}.dim(1),prior{im}.dim(3));
                    sippi_verbose(txt);
                end
            end
            
        end
    end
    
    %% SNESIM_STD OPTIONS (SGeMS version)
    if (strcmp(upper(prior{im}.type),'SNESIM_STD'))
        if ~isfield(prior{im},'S');
            prior{im}.S=sgems_get_par('snesim_std');
            sippi_verbose(sprintf('%s : Setting default SNESIM structure in prior{%d}.S',mfilename,im));
        else
            S=prior{im}.S;
            prior{im}.S=sgems_get_par('snesim_std');
            fn=fieldnames(S.XML.parameters);
            for ifn=1:length(fn);
                if isfield(S.XML.parameters,fn{ifn})
                    prior{im}.S.XML.parameters.(fn{ifn})=S.XML.parameters.(fn{ifn});
                end
            end           
            
            sippi_verbose(sprintf('%s : Updating with default SNESIM structure in prior{%d}.S',mfilename,im));
        end        
        prior{im}.S.dim.x=prior{im}.x;
        prior{im}.S.dim.y=prior{im}.y;
        prior{im}.S.dim.z=prior{im}.z;
        
        if isfield(prior{im},'ti');
            if isnumeric(prior{im}.ti);
                prior{im}.S.ti=prior{im}.ti;
            else
                prior{im}.S.ti_file=prior{im}.ti;
            end
            
        else
            prior{im}.ti=prior{im}.S.ti_file;
            sippi_verbose(sprintf('%s : Setting default training image as prior{%d}.ti=%s',mfilename,im,prior{im}.ti));
        end
    end
    
    %% SNESIM OPTIONS (V10 FORTRAN VERSION)
    if (strcmp(upper(prior{im}.type),'SNESIM'))
        if ~isfield(prior{im},'ti');
          prior{im}.ti=channels;
          sippi_verbose(sprintf('%s : Setting default training image as ''channels''',mfilename),0);
        end
        if ~isfield(prior{im},'S');
            prior{im}.S=snesim_init(prior{im}.ti,prior{im}.x,prior{im}.y,prior{im}.z);
            sippi_verbose(sprintf('%s : Setting default SNESIM structure in prior{%d}.S',mfilename,im),1);
        else
            prior{im}.S=snesim_init(prior{im}.ti,prior{im}.x,prior{im}.y,prior{im}.z,prior{im}.S);
            sippi_verbose(sprintf('%s : Updating SNESIM structure in prior{%d}.S',mfilename,im),2);
        end
        
        
        
        write_snesim(prior{im}.S);
        
    end
    
    
    
    
    %% SISIM OPTIONS
    if (strcmp(upper(prior{im}.type),'SISIM'))
        if ~isfield(prior{im},'S');
            prior{im}.S=sgems_get_par('sisim');
            sippi_verbose(sprintf('%s : Setting default SISIM structure in prior{%d}.S',mfilename,im));
        end
        prior{im}.S.dim.x=prior{im}.x;
        prior{im}.S.dim.y=prior{im}.y;
        prior{im}.S.dim.z=prior{im}.z;
            
        if ~isfield(prior{im},'Cm')&~isfield(prior{im},'Va');
            prior{im}.Cm='1 Sph(1)';
            txt=sprintf('%s : No covariance model set, using prior{%d}.Cm=''%s''',mfilename,im,prior{im}.Cm);
            sippi_verbose(txt);
        end       
        
        if ~isfield(prior{im},'marginal_prob')
            prior{im}.marginal_prob=prior{im}.S.XML.parameters.Marginal_Probabilities.value;
            txt_marginal_prob=sprintf('%g ',prior{im}.marginal_prob);
            txt=sprintf('%s : marginal probability not set. using prior{%d}.marginal_prob=[%s]',mfilename,im,txt_marginal_prob);
            sippi_verbose(txt);
        end
    end
    
    %% TARGET DIST - CHECK THAT NORMAL SCORE HASE BEEN PERFORMED
    if isfield(prior{im},'d_target')
        if ~iscell(prior{im}.d_target)
            if ~isfield(prior{im},'o_nscore')&&~iscell(prior{im}.d_target);
                [d_nscore,prior{im}.o_nscore]=nscore(prior{im}.d_target,1,1);
                sippi_verbose(sprintf('%s: performed normal score transformation of prior{%d}.d_target to prior{%d}.o_nscore',mfilename,im,im),1);
            end
        end
        if isfield(prior{im},'m0')
            if prior{im}.m0~=0
                sippi_verbose(sprintf('%s: Prior #%d, Both d_target and m0(>0, %g) have been set. Perhaps m0 should really be 0?',mfilename,im,prior{im}.m0))
                sippi_verbose(sprintf('%s: Prior #%d, -- This implies that the a prior mean is prior{%d}.m0+mean(prior{%d}.d_target)=%g',mfilename,im,im,im,prior{im}.m0+mean(prior{im}.d_target)))
                sippi_verbose(sprintf('%s: Prior #%d, Perhaps m0 should not be set, or be prior{%d}.m0=0; ',mfilename,im,im))
            end
        end
        
    end
   
    
    %%
    % OBSOLETE
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% WRITE TI IF APPLICABLE
    %if isfield(prior{im},'TI')
    %    sgems_write(prior{im}.S.ti_file,prior{im}.TI);
    %end
    
    %% confirmation of initialization;
    prior{im}.init=1;
    
    
end

