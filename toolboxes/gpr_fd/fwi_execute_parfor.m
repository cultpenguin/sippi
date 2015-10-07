function [dt nt error]=fwi_execute(Ncores_applied,ant_pos,sim_mode,addpar)

error=0;
tic
Positions=ant_pos;
Ncores=addpar.cores;

% Logicals which indicate if other then default values are applied in the
% simulation. Snapshot is always different from the default of the
% executable.
fr=0;
sn=1;
co=0; % 0=Default=Cartecial 
s1=2;

% Find out if an non-default value has been choosen
try
    coord=addpar.coord;
    if coord==1
        if addpar.debug==1
        disp('FDTD coordinate system: Cylindrical')
        end
        co=1;
    else
        if addpar.debug==1
        disp('FDTD coordinate system: Cartesian')
        end
    end
    
catch
    if addpar.debug==1
    disp('FDTD coordinate system: Cartesian (default)')
    end
end

try
    snap=addpar.snap;
    if addpar.debug==1
    disp(sprintf('%s %i','Snapshot "freqency" (in number of time-steps):',snap))
    end
catch
    snap=10000;
    if addpar.debug==1
    disp(sprintf('%s %i %s','Snapshot "freqency" (in number of time-steps):',10000,'~no snapshots (default)'))
    end
end

try
    frac=addpar.frac;
    fr=1;
    if addpar.debug==1
    disp(sprintf('%s %1.1f','Snapshot ratio:',frac))
    end
catch
    if addpar.debug==1
    disp(sprintf('%s %1.1f %s','Snapshot ratio:',1,'(default)'))
    end
end

% Determine the number of transmitter positions and where the associated
% receiver positions are located en the position array
try
    Ntrn=1;
    pos_new_trn(Ntrn)=1;
    for i=2:length(ant_pos(:,1))
        if ant_pos(i,2)~=ant_pos(i-1,2) || ant_pos(i,1)~=ant_pos(i-1,1)
            Ntrn=Ntrn+1;
            pos_new_trn(Ntrn)=i;
        end
    end
catch
    error=1;
    return
end

% Calculate the receiver positions used at the inverse model parameters:
%invfw_ratio=addpar.dx_inv/addpar.dx;
%z_pos=[(invfw_ratio/2)*addpar.dx:invfw_ratio*addpar.dx:addpar.Epsz*addpar.dx];
%x_pos=[(invfw_ratio/2)*addpar.dx:invfw_ratio*addpar.dx:addpar.Epsx*addpar.dx];

if sim_mode~=1    
    [X Z]=meshgrid(addpar.ipos_x*addpar.dx,addpar.ipos_z*addpar.dx);
    add_rec=[X(:) Z(:)];
    NREC=length(X(:));
end

% Start of the forward (and optional backward) simulation:
c_N_trn=0;
c_forward_data=0;
c_backward_data=0;
s2=1;
for i=1:ceil(Ntrn/Ncores_applied)
    c_N_cores=0;

    % Execute forward simulation
    for k=1:Ncores_applied
        c_N_cores=c_N_cores+1;
        c_N_trn=c_N_trn+1;
        try
            if c_N_trn<Ntrn
                trn=Positions(pos_new_trn(c_N_trn),1:2);
                rec=Positions(pos_new_trn(c_N_trn):pos_new_trn(c_N_trn+1)-1,3:4);
            else
                trn=Positions(pos_new_trn(c_N_trn),1:2);
                rec=Positions(pos_new_trn(c_N_trn):end,3:4);
                Ncores_applied=c_N_cores;
            end
        catch
            continue
        end
        
        % Write processing status to screen:
        if addpar.status==1
            disp(sprintf('Working on transmitter %d out of %d',c_N_trn,Ntrn));
            %progressbar(c_N_trn,Ntrn);
        end
        
        if sim_mode==2 || sim_mode==3
            rec=[rec;add_rec];
        end

        [fpath fname]=fileparts(which(eval(sprintf('%s%i','addpar.forwardexe',k))));
        fpath_names(k)=struct(['path' num2str(k)],{fpath});
        
        % Writing the trnrec.cor file containing all receiver positions
        delete(sprintf('%s%s%s',fpath,'\','*.bat'))
        delete(sprintf('%s%s%s',fpath,'\','*.dat'))
        delete(sprintf('%s%s%s',fpath,'\','*.cor'))
        
        ntrn=1;
        Fid=fopen(sprintf('%s%s',fpath,'\trnrec.cor'),'wt');
        fprintf ( Fid,'--------------------- Transmitter-Coordinates -----------------------\n');
        fprintf ( Fid,'%i\n',ntrn);
        for ii=1:ntrn
            fprintf(Fid,'%6.3f %6.3f\n',trn(ii,1),trn(ii,2));
        end;
        nrec(k)=length(rec(:,1));
        fprintf(Fid,'----------------------- Receiver-Coordinates ------------------------\n');
        fprintf(Fid,'%i\n',nrec(k));
        for ii=1:nrec(k)
            fprintf(Fid,'%6.3f %6.3f\n',rec(ii,1),rec(ii,2));
        end;
        fclose(Fid);

        % Setting up the batch-file
        fid=fopen(sprintf('%s%s',fpath,'\forward.bat'),'w');
        fprintf(fid,'%s\n','@echo off');
        fprintf(fid,'%s %s\n','cd',fpath);
        if co==0 && sn==0 && fr==0
            fprintf(fid,'%s%s\n',fname,'.exe');
        end
        if co==1 && sn==0 && fr==0
            fprintf(fid,'%s%s %s%i\n',fname,'.exe','/c:',coord);
        end
        if co==1 && sn==1 && fr==0
            fprintf(fid,'%s%s %s%i %s%i\n',fname,'.exe','/c:',coord,'/s:',snap);
        end
        if co==1 && sn==1 && fr==1
            fprintf(fid,'%s%s %s%i %s%i %s%i\n',fname,'.exe','/c:',coord,'/s:',snap,'/f:',frac);
        end
        if co==0 && sn==1 && fr==0
            fprintf(fid,'%s%s %s%i\n',fname,'.exe','/s:',snap);
        end
        if co==0 && sn==1 && fr==1
            fprintf(fid,'%s%s %s%i %s%i\n',fname,'.exe','/s:',snap,'/f:',frac);
        end
        if k<Ncores
            fprintf(fid,'%s','exit');
        end
        fclose(fid);
    end
                
    parfor (exe=1:Ncores_applied)   % Execute the forward FDTD simulation through forward.bat
        fpath=eval(['fpath_names.path' exe]);
        [s result]=dos(sprintf('%s%s',fpath,'\forward.bat'));
        
    end
    %if s==-1 
    %   disp('The external executable did not start')
    %        error=1;
    %        return
    %    end
    
    % Load the traces from the trace-files and move the content of interest
    % to the a mat-file in the working directory
    for j=1:Ncores_applied
        c_forward_data=c_forward_data+1;
        [fpath]=fileparts(which(eval(sprintf('%s%i','addpar.forwardexe',j))));
        
        % Load the output from the simulation
        % The program trys to load the product from the first core and so
        % fort. If e.g. core2 finishes before core1 the program has to wait
        % until core1 has finished.
        wait=1;
        while wait==1
            try
                Fid=fopen(sprintf('%s%s',fpath,'\trace_out000.dat'),'rb');
                tmp_nrec = fread(Fid,1,'int');
                nt = fread(Fid,1,'int');
                Trace_Er = fread(Fid,[nt,tmp_nrec],'double');
                Trace_Ez = fread(Fid,[nt,tmp_nrec],'double');
                Trace_Hphi = fread(Fid,[nt,tmp_nrec],'double');
                dt = fread(Fid,1,'double');
                fclose(Fid);
                wait=0;
            catch
                pause(0.001)
            end
        end


        clear Trace_Hphi

        %Trace=Trace_Ez((nrec(j)-NREC+1):end,:);
        if sim_mode==1 % Pure forward simulation
            switch addpar.output
                case{1}
                    clear Trace_Er
                    Trace_Ez=Trace_Ez';
                    save(sprintf('%s%i','observed_Ez_trn',c_forward_data),'Trace_Ez');
                case{2} 
                    clear Trace_Ez
                    Trace_Er=Trace_Er';
                    save(sprintf('%s%i','observed_Er_trn',c_forward_data),'Trace_Er');
                case{3} % 
                    Trace_Er=Trace_Er';
                    save(sprintf('%s%i','observed_Er_trn',c_forward_data),'Trace_Er');
                    Trace_Ez=Trace_Ez';
                    save(sprintf('%s%i','observed_Ez_trn',c_forward_data),'Trace_Ez');
                case{4} % As case=1 but the number on the output file is given by the vector addapr.no
                    if length(addpar.no)~=Ntrn 
                        disp('The length of addpar.no is different from Ntrn')
                        error=1;
                        return
                    end
                    clear Trace_Er
                    Trace_Ez=Trace_Ez';
                    save(sprintf('%s%i','observed_Ez_trn',addpar.no(c_forward_data)),'Trace_Ez');
                otherwise
                    disp('Warning: Invalid choise of addpar.output')
                    error=1;
                    return
            end
        elseif sim_mode==2 % Forward simulation as part of the fwi
            clear Trace_Er
            Trace_Ez=Trace_Ez';
            Trace_Ez_predict=Trace_Ez(1:nrec(j)-NREC,:);
            Trace_Ez=Trace_Ez((1+nrec(j)-NREC):end,:);
            save(sprintf('%s%i','forward_Ez_trn',c_forward_data),'Trace_Ez');
            save(sprintf('%s%i','predicted_Ez_trn',c_forward_data),'Trace_Ez_predict');
        elseif sim_mode==3 % Forward simulation for the perturbed model - part of the fwi 
            clear Trace_Er
            Trace_Ez=Trace_Ez';
            save(sprintf('%s%i','perturbed_Ez_trn',c_forward_data),'Trace_Ez');
        end

        if sim_mode==2
        
            d_obs=load(sprintf('%s%d%s','Observed data\observed_Ez_trn',c_forward_data,'.mat'),'-mat');

            % Unnormalised residuals:
            Residual_Ez=d_obs.Trace_Ez-Trace_Ez_predict;
            %Residual_Ez=Trace_Ez_resid;
            %Normalisation used during the calculation of the residuals for eps inversion:
%            X=[dt:dt:dt*nt];
%            tmp=Trace_Ez_resid;
%            faktor=sqrt(trapz(X,(d_obs.Trace_Ez').^2))./sqrt(trapz(X,(tmp').^2));
%            Residual_Ez=d_obs.Trace_Ez-(Trace_Ez_resid'*diag(faktor))';

            % Save the residuals in a binary - are afterwards loaded by backward.exe and backpropagated in time:
            [nx nz]=size(Residual_Ez);
            save2Dbin(sprintf('%s%s',fpath,'\resid00.E'),nx,nz,Residual_Ez');

            % The receiverpositions used in the backpropagation is stored
            % in a binary and read by backward.exe. Only receivers in the
            % inverse cells are nesecerry.
            Fid=fopen(sprintf('%s%s',fpath,'\resid00.N'),'wb');
            fprintf(Fid,'%d',nrec(j)-NREC);
            fclose(Fid);

        end
        
        % If snapshots are produced they are moved to the working directory
        if snap*dt<addpar.t
            for m=snap:snap:nt
                movefile(sprintf('%s\\EzS%05d_T00.dat',fpath,m),sprintf('forward_EzS%05d_T%02d.dat',m,c_forward_data-1))
            end
        end
    end

    if sim_mode==2
        delete(sprintf('%s%s%s',fpath,'\','*.dat'))
        % Execute backward simulation
        for k=1:Ncores_applied

            fpath=fileparts(which(eval(sprintf('%s%i','addpar.forwardexe',k))));
            fname='backward';

            % Setting up the batch-file
            fid=fopen(sprintf('%s%s',fpath,'\backward.bat'),'w');
            fprintf(fid,'%s\n','@echo off');
            fprintf(fid,'%s %s\n','cd',fpath);
            if co==0 && sn==0 && fr==0
                fprintf(fid,'%s%s\n',fname,'.exe');
            end
            if co==1 && sn==0 && fr==0
                fprintf(fid,'%s%s %s%i\n',fname,'.exe','/c:',coord);
            end
            if co==1 && sn==1 && fr==0
                fprintf(fid,'%s%s %s%i %s%i\n',fname,'.exe','/c:',coord,'/s:',snap);
            end
            if co==1 && sn==1 && fr==1
                fprintf(fid,'%s%s %s%i %s%i %s%i\n',fname,'.exe','/c:',coord,'/s:',snap,'/f:',frac);
            end
            if co==0 && sn==1 && fr==0
                fprintf(fid,'%s%s %s%i\n',fname,'.exe','/s:',snap);
            end
            if co==0 && sn==1 && fr==1
                fprintf(fid,'%s%s %s%i %s%i\n',fname,'.exe','/s:',snap,'/f:',frac);
            end
            if k<Ncores
                fprintf(fid,'%s','exit');
            end
            fclose(fid);

            % Execute the backward FDTD simulation through backward.bat
            if k<Ncores_applied
                s1=dos(sprintf('%s%s %s',fpath,'\backward.bat','&'));
            else
                s2=dos(sprintf('%s%s',fpath,'\backward.bat'));
            end

            if s1==-1 || s2==-1
                disp('The external backward executable did not start')
                error=1;
                return
            end
        end

        for j=1:Ncores_applied
            c_backward_data=c_backward_data+1;
            fpath=fileparts(which(eval(sprintf('%s%i','addpar.forwardexe',j))));

            % Load the output from the simulation
            % Load the output from the simulation
            % The program trys to load the product from the first core and so
            % fort. If e.g. core2 finishes before core1 the program has to wait
            % until core1 has finished.
            wait=1;
            while wait==1
                try
                    Fid = fopen(sprintf('%s%s',fpath,'\trace_out000.dat'),'rb');
                    tmp_nrec = fread(Fid,1,'int');
                    nt = fread(Fid,1,'int');
                    Trace_Er = fread(Fid,[nt,tmp_nrec],'double');
                    Trace_Ez = fread(Fid,[nt,tmp_nrec],'double');
                    Trace_Hphi = fread(Fid,[nt,tmp_nrec],'double');
                    dt = fread(Fid,1,'double');
                    fclose(Fid);
                    wait=0;
                catch
                    pause(0.5)
                end
            end

            clear Trace_Er
            clear Trace_Hphi
            Trace_Ez=Trace_Ez';

            Trace_Ez=Trace_Ez(nrec(j)-NREC+1:end,:);
            save(sprintf('%s%i','backward_Ez_trn',c_backward_data),'Trace_Ez');

            if snap*dt<addpar.t
                for m=snap:snap:nt
                    movefile(sprintf('%s\\EzS%05d_T00.dat',fpath,m),sprintf('backward_EzS%05d_T%02d.dat',m,c_backward_data-1))
                end
            end
        end
    end
end

time=toc;
if time<1
    disp('The external executable did not start')
    error=1;
    return
end

hr=floor(time/3600);
min=floor((time/60-hr*60));
sec=floor((time-min*60-hr*3600));

if addpar.debug==1
disp('-----------------------------------------------')
disp(sprintf('%s %i %s %i %s %i %s','The FDTD scheme terminated successfully after:',hr,'hr,',min,'min, and',sec','sec'))
disp('-----------------------------------------------')
end

function save2Dbin(fname,nx,nz,data)
Fid=fopen(sprintf('%s',fname),'wb');
fwrite(Fid,nx,'int');
fwrite(Fid,nz,'int');
fwrite(Fid,data,'double');
fclose(Fid);