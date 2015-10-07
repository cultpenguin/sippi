%==========================================================================
%The function [TRC,CFG] = LOADTRACEM(FIELD,TRC_REDUCTION,SRC_NR,
%                                    TOOLTYPE) 
%loads traces from trace-files, created with the FDTD algorithm fwi.
%
%FIELD         = field that should be loaded ('-1'=all -> returns a
%                structure or single field components '1'=Er, '2'=Ez or
%                '3'=Hphi(NEEDED)
%TRC_REDUCTION = trace-length 'nt' (#samples) to load. If set to '-1',
%                everything is loaded.
%                CAUTION: if TOOLTYPE==4x, this denotes the iteration 
%                number (NEEDED)
%SRC_NR        = shot-number to load (NEEDED)
%TOOLTYPE      = if 1=cartesian modeling tool; 2=cylindrical modeling tool,
%                3=reflex-data and 4x=cart+wave forminversion, where x
%                stands for 0 (simulated-), 1 (perturbed-), 2 (residual-) and 
%                3 (observed-) traces. (NEEDED)
%ECHOTYPE      = either 1 (every comment is ploted or 0 to disable the
%                comments). (NOT NEEDED)
%
%The return parameters are: Trc containing the trace(s); and cfg containing
%the configuration parameters
%
%REMARK: Concerning the Waveform-Inversion - the residual traces after an
%Epsilon-Inversion are ALLWAYS scaled in order to get rid of the Sigma-
%effects! So if REAL residual traces are needed, compute them by simply
%subtracting the observed from the simulated traces!
%
%(c) Jacques Ernst                                               ETH Zurich
%==========================================================================
function [Trc,cfg] = LoadTraceM(field,trc_reduction,src_nr,tooltype,echotype)


echotype=-1;

if (tooltype>=40) 
    fieldtype=tooltype-40;
    tooltype=4;
end;


%Load Config DATA
switch tooltype
    case {2} %Cyli. Coordinate Tool
        Fid = fopen(sprintf('%s','radarin.dat'),'rb');
        if (Fid==-1),
            uiwait(errordlg('Config-File not found','File Error','modal'));
            Trc=[];cfg=[];
            return
        else
            cfg.ver = fread(Fid,1,'int');
            cfg.dr = fread(Fid,1,'double');
            fread(Fid,1,'int'); %Nr_h
            if (cfg.ver < 8) fread(Fid,1,'int'); end;%zero
            fread(Fid,1,'int'); %Nr_l
            fread(Fid,1,'int'); %Nz_l
            cfg.border = fread(Fid,1,'int');
            if (cfg.ver > 5), fread(Fid,1,'double'); end;
            cfg.gather = fread(Fid,1,'int');
            cfg.fc = fread(Fid,1,'double');
            fread(Fid,1,'double'); %T
            cfg.dt = fread(Fid,1,'double'); %Dt
            cfg.vmax = fread(Fid,1,'double');
            cfg.nsrc = fread(Fid,1,'int');
            cfg.Sr = fread(Fid,[1,cfg.nsrc],'double');
            cfg.Sz = fread(Fid,[1,cfg.nsrc],'double');
            cfg.nrec = fread(Fid,1,'int');
            cfg.Rr = fread(Fid,[1,cfg.nrec],'double');
            cfg.Rz = fread(Fid,[1,cfg.nrec],'double');
            fread(Fid,1,'int'); %wavelet
            if (cfg.ver > 7) fread(Fid,1,'int'); end; %fixed/variable source
            fread(Fid,1,'int'); %snap
            fread(Fid,1,'int'); %fout
            cfg.a_type = fread(Fid,1,'int'); %antenna type
            cfg.a_r = fread(Fid,1,'double'); %antenna radius
            cfg.a_l = fread(Fid,1,'double');
            cfg.a_ddz = fread(Fid,1,'double'); %antenna gap
            fclose(Fid);
        end;
end;
%Load Trace DATA
%Load Config DATA
switch tooltype 
    case {1,2} %Cart. or Cyli. Coordinate Tool
        if ~exist(sprintf('%s','trace_out00.dat')),
            Fid = fopen(sprintf('trace_out%03d.dat',src_nr),'rb');
            if echotype,
                %disp(sprintf('\nLoading %s\\trace_out%03d.dat\n',path,src_nr));
            end;
        else
            Fid = fopen(sprintf('%s','trace_out00.dat'),'rb');
            if echotype,
                %disp(sprintf('\nLoading %s\\trace_out00.dat\n',path));
            end;
        end;
        if (Fid==-1),
            uiwait(errordlg('Trace-File not found','File Error','modal'));
            Trc=[];cfg=[];
            return
        else
            cfg.nrec = fread(Fid,1,'int');
            cfg.nt = fread(Fid,1,'int');
            Trace_Er = fread(Fid,[cfg.nt,cfg.nrec],'double');
            Trace_Ez = fread(Fid,[cfg.nt,cfg.nrec],'double');
            Trace_Hphi = fread(Fid,[cfg.nt,cfg.nrec],'double');
            if (tooltype==1), cfg.dt = fread(Fid,1,'double'); end;
            fclose(Fid);

            if trc_reduction ~= -1,
                nt_temp = cfg.nt;
                cfg.nt = trc_reduction;
                if echotype,
                    %disp(sprintf('\nWARNING: Trace reduction applied!! Trace length is now: %i instead of %i samples!!\n',...
                     %   trc_reduction,nt_temp));
                end;
            end;
        end;
    case {3} %REFLEX
        D = dir(sprintf('%s','*__.RAD'));
        if isempty(D), Trc=[];cfg=[]; uiwait(errordlg('Trace-File not found','File Error','modal')); return; end;
        if echotype,
            %disp(sprintf('\nLoading %s corresponding to shot number %i\n',...
            %    D(src_nr).name,src_nr));
        end;

        cfg.ver = 99999; cfg.dr = 99999; cfg.border = 99999; 
        cfg.gather = 99999; cfg.fc = 99999; cfg.vmax = 99999;
        cfg.ntrn = length(D);
        cfg.Sr = 99999; cfg.Sz = 99999;

        warning off MATLAB:divideByZero
        [cfg.nt,cfg.nrec,cfg.dt,Trace_Ez,dx_rec] = ...
            ReadOneProf(sprintf('\%s',D(src_nr).name(1:find(D(src_nr).name=='.')-1)));
        cfg.dt = cfg.dt*1e-9; %ns -> s
        field = 2; % only Ez fields are saved!

        cfg.Rr = 99999;
        cfg.Rz = (1:cfg.nrec).*dx_rec;
        cfg.a_type = 99999; cfg.a_r = 99999; cfg.a_l = 99999;
        cfg.a_ddz = 99999;
    case {4} %cart. coordinate using fwinversions code!!
        Fid = fopen(sprintf('\traces_IT%02d_D%01i.dat',trc_reduction,fieldtype),'rb');
        if echotype,
            %disp(sprintf('%s\\traces_IT%02d_D%01i.dat',trc_reduction,fieldtype));
        end;

        if (Fid==-1),
            uiwait(errordlg('Trace-File not found','File Error','modal'));
            Trc=[];cfg=[];
            return
        else
            trc_type = fread(Fid,1,'int');
            cfg.ntrn = fread(Fid,1,'int');
            if (src_nr>=cfg.ntrn)
                disp('WARNING: No valid transmitter number! Exiting!!');
                Trc=[];cfg=[];
                return;
            end;
            cfg.nrec = fread(Fid,1,'int');
            cfg.dt = fread(Fid,1,'double');
            cfg.nt = fread(Fid,1,'int');
            fseek(Fid,(src_nr*cfg.nrec*cfg.nt)*8,'cof');               
            tmp_E = fread(Fid,cfg.nt*cfg.nrec,'double');
            %Eobs = fread(Fid,ntrn*cfg.nt*cfg.nrec,'double');
            %fseek(Fid,-8,'eof');
            fclose(Fid);

            field = 2; % only Ez fields are saved!

            Trace_Ez = reshape(tmp_E,cfg.nt,cfg.nrec);

            clear tmp_E;
            
            if echotype,
                if trc_reduction ~= -1,
                    %disp(sprintf('\nINFO: Traces for iteration %i are shown!!\n',...
                     %   trc_reduction));
                end;

                switch (trc_type)
                    case {0}
                        outp = 'Synthetic Data';
                    case {1}
                        outp = 'Perturbed Data';
                    case {2}
                        outp = 'Data Residuals';
                    case {3}
                        outp = 'Observed Data';
                end;
                disp(sprintf('Output is: %s\n',outp));
            end;
        end;
    case {5} %Sensitivite Trace-Files (not Sensitivites!!)
        if (field==1) %Er
            fieldtype = 1;
        elseif (field==2) %Ez
            fieldtype = 0;
        end;
        Fid = fopen(sprintf('\Sens%03d_F%01i.dat',src_nr,fieldtype),'rb');
        if echotype,
            disp(sprintf('\Sens%03d_F%01i.dat',src_nr,fieldtype));
        end;

        if (Fid==-1),
            uiwait(errordlg('Trace-File not found','File Error','modal'));
            Trc=[];cfg=[];
            return
        else
            cfg.nrec = fread(Fid,1,'int');
            cfg.nt   = fread(Fid,1,'int');
            Trace_Ez = fread(Fid,[cfg.nt,cfg.nrec],'double');
            cfg.dt   = fread(Fid,1,'double');
            fclose(Fid);
        end;
end;

switch field
    case {-1}
        if echotype,
            disp('All field components are loaded!');
        end;
        Trc.Er = Trace_Er(1:cfg.nt,:);
        Trc.Ez = Trace_Ez(1:cfg.nt,:);
        Trc.Hphi = Trace_Hphi(1:cfg.nt,:);
        cfg.M_Er = max(abs(Trc.Er));
        cfg.m_Er = min(abs(Trc.Er));
        cfg.M_Ez = max(abs(Trc.Ez));
        cfg.m_Ez = min(abs(Trc.Ez));
        cfg.M_Hphi = max(abs(Trc.Hphi));
        cfg.m_Hphi = min(abs(Trc.Hphi));
    case {1}
        if echotype,
            %disp('Only the Ex/Er-field component is loaded!');
        end;
        Trc = Trace_Er(1:cfg.nt,:);
        cfg.M_Er = max(abs(Trc));
        cfg.m_Er = min(abs(Trc));
        cfg.M_Ez = 99999;
        cfg.m_Ez = 99999;
        cfg.M_Hphi = 99999;
        cfg.m_Hphi = 99999;
    case {2}
        if echotype,
            disp('Only the Ez-field component is loaded!');
        end;
        Trc = Trace_Ez(1:cfg.nt,:);
        cfg.M_Er = 99999;
        cfg.m_Er = 99999;
        cfg.M_Ez = max(abs(Trc));
        cfg.m_Ez = min(abs(Trc));
        cfg.M_Hphi = 99999;
        cfg.m_Hphi = 99999;
    case {3}
        if echotype,
            disp('Only the Hy/Hphi-field component is loaded!');
        end;
        Trc = Trace_Hphi(1:cfg.nt,:);
        cfg.M_Er = 99999;
        cfg.m_Er = 99999;
        cfg.M_Ez = 99999;
        cfg.m_Ez = 99999;
        cfg.M_Hphi = max(abs(Trc));
        cfg.m_Hphi = max(abs(Trc));
end;
clear Trace_Er Trace_Ez Trace_Hphi;
