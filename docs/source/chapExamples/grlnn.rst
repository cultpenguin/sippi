Cross hole GPR tomopgraphy using Neural Networks
================================================

The example below requires the `SIPPI
toolbox <http//sippi.sourceforge.net/>`__.

::

    % grl_nn 
    clear all;close all

    createTrainingSet=1;
    createReferenceModel=1;
    useTargetDist=0;
    Ntrain=40000;

    Nr_modeling=6000; %size of sample for modeling error

    TrainSizes=[1000 5000 10000 20000 40000];

    splitData=3;
    epochs=30000; 
    hiddenLayerSize=80;

    %
    %Ntrain=2500;
    %splitData=0;
    %epochs=10000; 
    %hiddenLayerSize=20;
    %TrainSizes=[1000 2500];


    % mcmc 
    nite=2000000;;

    %% LOAD DATA AND CONFIGURATION
    D=load('AM13_data.mat');

    % REVERSE S and R for fd forward
    D2=D;
    D.S(352:end,:)=D2.R(352:end,:);
    D.R(352:end,:)=D2.S(352:end,:);
    D.S(:,1)=D.S(:,1)+1;
    D.R(:,1)=D.R(:,1)+1;
    clear D2
    i_use=1:1:size(D.S,1);
    id=1;
    data{id}.d_obs=D.d_obs(i_use);
    data{id}.d_std=D.d_std(i_use).*0+0.1;


    %% B: Define FD
    %forward_fd;
    forward.forward_function='sippi_forward_traveltime';
    forward.sources=D.S(i_use,:);
    forward.receivers=D.R(i_use,:);
    forward.type='fd';
    %forward.type='fat';forward.linear=1;forward.freq=0.1;
    %forward.type='eikonal';
    %forward.m0=prior{1}.m0;

    forward.fd.t=1e-7;
    forward.fd.addpar.Tg=100*10^6;
    forward.fd.dx_fwd=0.1;

    %% A: Define prior model
    im=1;
    dx=0.2;
    prior{im}.type='FFTMA';
    prior{im}.name='Velocity (m/ns)';
    prior{im}.m0=0.1431;
    prior{im}.Va='.000215 Sph(6)';
    prior{im}.x=[-.5:dx:6.5];
    prior{im}.x=[0:dx:7];
    prior{im}.y=[0:dx:13];
    prior{im}.cax=[-1 1].*.03+prior{im}.m0;

    if useTargetDist==1;
        d_target=[randn(1,100)*.003-0.01 randn(1,100)*.003+0.01]+prior{im}.m0;
        prior{im}.d_target=d_target;
        prior{im}.m0=0; %% MAKE SURE sippi_forward_traveltime tests for a non-zero velocity
    end


    %% C: Make Reference model
    if createReferenceModel==1
        rng(1);
        % Reference mo
        [m_ref,prior]=sippi_prior(prior);
        NM=prod(size(m_ref{1}));
        % reference data
        [d_ref,forward]=sippi_forward(m_ref,forward,prior);

        % data    
        data{1}.d_ref=d_ref{1};
        data{1}.d_noise=randn(size(d_ref{1})).*data{1}.d_std;
        data{1}.d_obs=data{1}.d_ref+data{1}.d_noise;



        % compute reference data in m0
        m_ref0=m_ref;
        m_ref0{1}=m_ref0{1}.*0+prior{1}.m0;
        disp('computing reference forward')
        [d_ref0,forward0]=sippi_forward(m_ref0,forward,prior);
        d0=d_ref0{1};

        save grl_ReferenceModel 
    else
        load grl_ReferenceModel
    end

    %% D: Create M-D training data set for forward model
    if createTrainingSet==1
        ATTS=zeros(length(m_ref{1}(:)),Ntrain);
        DATA=zeros(length(d_ref{1}(:)),Ntrain);

        iplot=1;
        t0=now;
        for i=1:Ntrain;

            if ((i/iplot)==round(i/iplot)&&(i>1));progress_txt(i,Ntrain,time_loop_end(t0,i-1,Ntrain));end
            m=sippi_prior(prior);
            try
                d=sippi_forward(m,forward,prior,data);

                ATTS(:,i)=m{1}(:);
                DATA(:,i)=d{1}(:)-d0;
            catch
                disp(sprintf('Something went wrong.'));
                i=i-1;
            end
        end
        save(sprintf('grl_%s_NM%d_NT%d',forward.type,NM,Ntrain));
    else
        load grl_eikonal_NM2376_NT300
        %load grl_eikonal_NM2376_NT1000.mat
    end


    %% E: SETUP FORWARD MODELS

    if ~exist('splitData'); 
        splitData=0; % SPLIT DATA FOR NN
    end
    if ~exist('epochs'); 
        epochs=100000; %   
    end
    if ~exist('hiddenLayerSize'); 
        hiddenLayerSize=80;
    end

    forward_nn.forward_function='sippi_forward_mynn';
    forward_nn.sources=forward.sources;
    forward_nn.receivers=forward.receivers;
    forward_nn.ATTS=ATTS;
    forward_nn.DATA=DATA;
    forward_nn.d0=d0;
    clear DATA ATTS
    forward_nn.splitData=splitData;
    forward_nn.epochs=epochs;;  
    forward_nn.hiddenLayerSize=hiddenLayerSize;
    forward_nn.max_nm=1e+10;
    txt=sprintf('grl_NM%d_DX%d_%s_NT%d_SD%d_NH%d',NM,dx*100,forward.type,Ntrain,forward_nn.splitData,forward_nn.hiddenLayerSize);
    forward_nn.mfunc_string=txt;
    disp(sprintf('%s: setting name ''%s''',mfilename,txt));

    % setup all forward models
    i_forward=0;
    if ~exist('TrainSizes'); 
        TrainSizes=[100 200];
    end
    for n_use=TrainSizes;
        if n_use<= size(forward_nn.ATTS,2);
            i_forward=i_forward+1;
            f_mul{i_forward}=forward_nn;
            f_mul{i_forward}.ATTS=forward_nn.ATTS(:,1:n_use);
            f_mul{i_forward}.DATA=forward_nn.DATA(:,1:n_use);       
            txt_use=sprintf('grl_NM%d_DX%d_%s_NT%d_SD%d_NH%d',NM,dx*100,forward.type,n_use,forward_nn.splitData,forward_nn.hiddenLayerSize);
            f_mul{i_forward}.mfunc_string=txt_use;
        end
    end

    % eikonal
    i_forward=i_forward+1;
    f_mul{i_forward}=forward;
    f_mul{i_forward}.type='eikonal';

    % ray_2d
    i_forward=i_forward+1;
    f_mul{i_forward}=forward;
    f_mul{i_forward}.type='ray_2d';
    f_mul{i_forward}.linear=1;
    f_mul{i_forward}.freq=0.1;
    f_mul{i_forward}.r=2;
    f_mul{i_forward}.normalize_vertical=0;


    %% F: EVALUATE forward models once to setup NN and Linear operators
    for i=1:length(f_mul);
        t1=now;
        [d_mul{i},f_mul{i}]=sippi_forward(m_ref,f_mul{i},prior);
        t2=now;
        time_mul{i}=(t2-t1)*3600*24;

    end

    save(sprintf('%s_forward',txt))

    %% G: Estimate modeling errors
    if ~exist('Nr_modeling'); 
        Nr_modeling=6000;
    end
    [Ct,dt,dd,d_full,d_app]=sippi_compute_modelization_forward_error(forward,f_mul,prior,Nr_modeling);

    % H: Setup one data structure per forward model, with the correct modeling error
    for i=1:length(f_mul);    

        %s=sum(abs(dd{5}));
        %ii=find(s<180);

        data_mul{i}=data;
        data_mul{i}{1}.dt=dt{i};
        data_mul{i}{1}.Ct=Ct{i};
    end

    save(sprintf('%s_modelerr',txt))


    %% I: Perform probabilistic inversion using extended Metropolis sampling
    if ~exist('nite'); 
        nite=1000000; % 
    end
    options.mcmc.m_ref=m_ref;
    options.mcmc.nite=nite;   % [1] : Number if iterations
    options.mcmc.i_sample=ceil(options.mcmc.nite/1000); % : Number of iterations between saving model to disk
    options.mcmc.i_plot=100000;  % [1]: Number of iterations between updating plots
    options.mcmc.i_save_workspace=100000;  % [1]: Number of iterations between

    i_burnin=options.mcmc.nite/30;
    prior{1}.seq_gibbs.i_update_step_max=i_burnin;

    options.mcmc.anneal.i_begin=1; % default, iteration number when annealing begins
    options.mcmc.anneal.i_end=ceil(i_burnin/2); %  iteration number when annealing stops
    options.mcmc.anneal.T_begin=5; % Start temperature for annealing
    options.mcmc.anneal.T_end=1; % End temperature for annealing

    %options.mcmc.n_chains=2; % set number of chains (def=1)
    %options.mcmc.T=[1 1.1 1.2];      % set temperature of chains [1:n_chains]

    % RUN MCMC

    rseed=1;
    for i=1:(length(f_mul));
        rng(rseed);
        if isfield(f_mul{i},'mfunc');
            options.txt=f_mul{i}.mfunc_string;        
        else
            try
                %options.txt=sprintf('grl_NM%d_DX%d_%s_SD%d_NH%d',NM,dx*100,f_mul{i}.type,forward_nn.splitData,forward_nn.hiddenLayerSize);
                options.txt=sprintf('grl_NM%d_DX%d_%s_SD%d_NH%d',NM,dx*100,forward_nn.type,forward_nn.splitData,forward_nn.hiddenLayerSize);
            catch
                options.txt=txt;
            end
        end
        t_start=now;
        [o]=sippi_metropolis(data_mul{i},prior,f_mul{i},options);
        options_out{i}.txt=o.txt;
        %sippi_plot_posterior_sample(options_out{i}.txt);
        sim_minutes(i)=(now-t_start)*60*24;
    end

    save(sprintf('%s_inverted',txt))
