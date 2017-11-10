% sippi_forward_mynn: Neural Network forward
%  
% A training data set consisting on N sets of model parameters and data 
% needs to available as [ATTS,DATA]
%
%
% % required fields
% forward.ATTS  [NM,N]
% forward.DATA [ND,N]
% where 
%   NM = length(m{1}(:))
%   ND = length(data{1}.d_obs(:))
%
% Optional:
%    forward.epochs=50000;       % number of iterations forback propagation
%    forward.hiddenLayerSize=10; % number of hidden layers
%
%
% See also: sippi_forward, neural_network_learn
function [d,forward,prior,data]=sippi_forward_mynn(m,forward,prior,data,id,im)


if ~isfield(forward,'xx');
    [forward.xx,forward.yy]=meshgrid(prior{1}.x,prior{1}.y);
end

if ~isfield(forward,'mfunc_string');
    forward.mfunc_string='nn';
end

if ~isfield(forward,'epochs');
    forward.epochs=50000;;
end
if ~isfield(forward,'hiddenLayerSize');
    forward.hiddenLayerSize=10;;
end

if size(forward.ATTS,2)<25000
    useGPU=1;
else
    useGPU=0;
end
trainParam.epochs=forward.epochs;
                

if isfield(forward,'mfunc');
    % solve forward problem
    if length(forward.mfunc)==1;
        try
            d{1}=feval(forward.mfunc{1},m{1}(forward.im{1}));
        catch
            d{1}=feval(forward.mfunc{1},m{1}(:));
        end
    else
        nd=size(forward.DATA,1);
        d{1}=NaN.*ones(nd,1);
        for i=1:length(forward.net)
            % evaluate using function 
            d{1}(forward.id{i})=feval(forward.mfunc{i},m{1}(forward.im{i}(:)));
            % evaluate using NET structure
            %d{1}(forward.id{i})=forward.net{i}(m{1}(forward.im{i}(:)));
        end
    end
else
    
    if isfield(forward,'net')
        % evaluate net
    else
        % setup net
        if ~isfield(forward,'splitData');
            forward.splitData=2;
        end
        
        if forward.splitData==1,
            
            j=0;
            id_rest=[1:1:size(forward.DATA,1)];
            for minY=[0,4,8]
                j=j+1;
                maxY=minY+4;
                forward.id{j}=find(forward.sources(:,2)<maxY & forward.receivers(:,2)<maxY & forward.sources(:,2)>=minY & forward.receivers(:,2)>=minY);
                forward.im{j}=find(forward.yy<maxY & forward.yy>minY);
                
                forward.mfunc{j}=sprintf('%s_%02d',forward.mfunc_string,j);
                [forward.net{j},tr,forward.mfunc{j}]=traveltime_nn_learn_neural_network(forward.ATTS(forward.im{j},:), forward.DATA(forward.id{j},:),forward.mfunc{j});
                id_rest=setxor(id_rest,forward.id{j});
                
            end
            
            j=j+1;
            forward.im{j}=1:prod(size(forward.xx));
            forward.im{j}=forward.im{j}(:);
            forward.id{j}=id_rest;
            
            forward.mfunc{j}=sprintf('%s_%02d',forward.mfunc_string,j);
            [forward.net{j},tr,forward.mfunc{j}]=traveltime_nn_learn_neural_network(forward.ATTS(forward.im{j},:), forward.DATA(forward.id{j},:),forward.mfunc{j},forward.hiddenLayerSize,trainParam,useGPU);
            
            
        elseif forward.splitData==2,
            % equal chunks
            %d_sr=40;
            d_sr=360;
            nd=size(forward.DATA,1);
            j=0;
            for i1=1:d_sr:nd;
                i2=min([i1+d_sr-1 nd]);
                disp(sprintf('%s: learning NN for data %d:%d',mfilename,i1,i2));
                j=j+1;
                
                % use subset of model parameters.. perhaps use wavelets?
                im{j}=1:prod(size(forward.xx));
                nm=length(im{j});
                max_m=360;
                im_use=ceil(linspace(1,nm,min([nm max_m])));
                
                im{j}=im{j}(im_use);
                id{j}=i1:1:i2;
                mfunc{j}=sprintf('%s_%02d',forward.mfunc_string,j);
                try
                    [forward.net{j},tr,mfunc{j}]=traveltime_nn_learn_neural_network(forward.ATTS(im{j},:), forward.DATA(id{j},:),mfunc{j},forward.hiddenLayerSize,trainParam,useGPU);
                catch
                    keyboard
                end
                
            end
            
            forward.mfunc=mfunc;
            forward.id=id;
            forward.im=im;
        elseif forward.splitData==3;
            %% SPLIT DATA IN CHUNBKS WITH SIMULAR SourcY adn ReceiverY locations, 
            %  AND USE A SUBSET OF MODEL PARAMETERS SPANNED BY THE S-R
            
            if ~isfield(forward,'max_nm');
                forward.max_nm=1e+9;
            end
            
            j=0;
            id_rest=[1:1:size(forward.DATA,1)];
            dY=1;
            ddY=max([0 5.5-dY]);
            for Y=[0:dY:12];
                minY=Y-ddY;
                maxY=(Y+dY+ddY);
                id=find(forward.sources(:,2)<maxY & forward.receivers(:,2)<maxY & forward.sources(:,2)>=minY & forward.receivers(:,2)>=minY);;
                id=intersect(id,id_rest);
                
                
                if length(id)>0
                    minSY = min([ min(forward.sources(id,2)) min(forward.receivers(id,2)) ])-.5;
                    maxSY = max([ max(forward.sources(id,2)) max(forward.receivers(id,2)) ])+.5;
                    im=find(forward.yy<=maxSY & forward.yy>=minSY);
                    
                    nm=length(im);
                    %max_nm=360;
                    %im_use=ceil(linspace(1,nm,min([nm max_nm])));
                    im_use=randomsample(1:nm,min([nm forward.max_nm]));
                    im=im(im_use);
                    
                    
                    j=j+1;
                    forward.id{j}=id;
                    forward.im{j}=im;
                    
                    id_rest=setxor(id_rest,forward.id{j});
                end
            end
            
            if ~isempty(id_rest)
                j=j+1;
                forward.id{j}=id_rest;
                forward.im{j}=1:prod(size(forward.yy));
            end
                
            %%
            for j=1:length(forward.id);
            
                figure(20);
                subplot(1,10,j)
                imagesc(prior{1}.x,prior{1}.y,m{1});
                set(gca,'FontSize',6);
                axis image
                hold on
                plot(forward.xx(forward.im{j}),forward.yy(forward.im{j}),'g.');
                plot([forward.sources(forward.id{j},1),forward.receivers(forward.id{j},1)]',[forward.sources(forward.id{j},2),forward.receivers(forward.id{j},2)]','k-*')
                hold off
                drawnow
                print_mul(forward.mfunc_string);
                
                forward.mfunc{j}=sprintf('%s_%02d',forward.mfunc_string,j);
                
                %[forward.net{j},tr,forward.mfunc{j}]=traveltime_nn_learn_neural_network(forward.ATTS(forward.im{j},:), forward.DATA(forward.id{j},:),forward.mfunc{j},forward.hiddenLayerSize,forward.epochs);                
                [forward.net{j},tr,forward.mfunc{j}]=neural_network_learn(forward.ATTS(forward.im{j},:), forward.DATA(forward.id{j},:),forward.mfunc{j},forward.hiddenLayerSize,trainParam,useGPU);                

            end
        elseif forward.splitData>10;
            % USE NO SUBSET OF DATA. ONLY SUBSET OF MODELS
            
            nd=size(forward.DATA,1);
            id{1}=1:1:nd;
            
            dm=ceil(forward.splitData-10);
            sippi_verbose(sprintf('%s: Using every ''%d'' pixel in every dimension.',mfilename,dm));
            ii=forward.xx.*0;
            ii(1:dm:end,1:dm:end)=1;
            im{1}=find(ii==1)';
            %im{1}=1:prod(size(forward.xx));
            im{1}=im{1}(:);
            
            
            forward.id{1}=id{1}(:);
            forward.im{1}=im{1}(:);
            
            forward.mfunc{1}=sprintf('%s_dm%d',forward.mfunc_string,dm);
            [forward.net{1},tr,mfunc{1}]=neural_network_learn(forward.ATTS(forward.im{1},:), forward.DATA(forward.id{1},:),forward.mfunc{1},forward.hiddenLayerSize,trainParam,useGPU);
            
            %forward.mfunc=mfunc;
            %forward.id=id;
            %forward.im=im;
            
            
        else
            % USE ALL DATA AND ALL MODELPARAMETERS
            nd=size(forward.DATA,1);
            id{1}=1:1:nd;
            im{1}=1:prod(size(forward.xx));
            im{1}=im{1}(:);
            
            mfunc{1}=sprintf('%s',forward.mfunc_string);
            [forward.net{1},tr,mfunc{1}]=neural_network_learn(forward.ATTS, forward.DATA,mfunc{1},forward.hiddenLayerSize,trainParam,useGPU);
            
            forward.mfunc=mfunc;
            forward.id=id;
            forward.im=im;
            
        end
        %DATA_SIM = forward.net(forward.ATTS);
        %forward.performance = perform(forward.net,forward.DATA,DATA_SIM);
        %forward.mfunc=mfunc;
        
    end
    
    % Solve forward problem    
    nd=size(forward.DATA,1);
    d{1}=NaN.*ones(nd,1);
    for i=1:length(forward.net),
        d{i}(forward.id{i}(:))=forward.net{i}(m{1}(forward.im{i}(:)));
    end
end

% network predicts variations around d0
if isfield(forward,'d0'),
    d{1}=d{1}+forward.d0;
end

