
function [addpar error]=write_parameters_to_screen(ant_pos,Sig,Eps,addpar)

error=0;

try
    %=============== Parameters used in the forward modelling ================%
    disp('-----------------------------------------------')

    % Name of forward controler #1:
    try
        disp(sprintf('%s %s','Name of forward executable #1:',addpar.forwardexe1))
    catch
        addpar.forwardexe1='FDTD_forward1.exe';
        disp(sprintf('%s %s %s','Name of forward executable #1:',addpar.forwardexe1,'(default)'))
    end

    % Name of forward controler #2:
    try
        disp(sprintf('%s %s','Name of forward executable #2:',addpar.forwardexe2))
    catch
        addpar.forwardexe2='FDTD_forward2.exe';
        disp(sprintf('%s %s %s','Name of forward executable #2:',addpar.forwardexe2,'(default)'))
    end

    % Name of forward controler #3:
    try
        disp(sprintf('%s %s','Name of forward executable #3:',addpar.forwardexe3))
    catch
        addpar.forwardexe3='FDTD_forward3.exe';
        disp(sprintf('%s %s %s','Name of forward executable #3:',addpar.forwardexe3,'(default)'))
    end

    % Name of forward controler #4:
    try
        disp(sprintf('%s %s','Name of forward executable #4:',addpar.forwardexe4))
    catch
        addpar.forwardexe4='FDTD_forward4.exe';
        disp(sprintf('%s %s %s','Name of forward executable #4:',addpar.forwardexe4,'(default)'))
    end

    disp('-----------------------------------------------')

    % Number of boundary cells:
    try
        disp(sprintf('%s %i','Number of boundary cells:',addpar.bx))
    catch
        addpar.bx=40;
        disp(sprintf('%s %i %s','Number of boundary cells:',addpar.bx,'(default)'))
    end

    % GPML-Tuning: sm (-1:auto): Max conductivity:
    try
        disp(sprintf('%s %f','GPML-Tuning (S_max):',addpar.sm))
    catch
        addpar.sm=-1;
        disp(sprintf('%s %s %s','GPML-Tuning (S_max):','auto','(default)'))
    end

    % GPML-Tuning: sigm (-1:auto): Max coord. strecking:
    try
        disp(sprintf('%s %f','GPML-Tuning (Sig_max):',addpar.sigm))
    catch
        addpar.sigm=-1;
        disp(sprintf('%s %s %s','GPML-Tuning (Sig_max):','auto','(default)'))
    end

    % GPML-Tuning expo: Exponent, describing how fast the strecting and
    % conductivity happens in the boundaries:
    try
        disp(sprintf('%s %i','GPML-Tuning (Exponent):',addpar.n))
    catch
        addpar.n=2;
        disp(sprintf('%s %i %s','GPML-Tuning (Exponent):',addpar.n,'(default)'))
    end

    disp('-----------------------------------------------')
    % C-Stability-Scaling (do not change): Stability modifier, controls
    % dispersion criteria (Courant):
    try
        disp(sprintf('%s %f','Courant stability scaling:',addpar.sc))
    catch
        addpar.sc=0.9;
        disp(sprintf('%s %f %s','Courant stability scaling:',addpar.sc,'(default)'))
    end

    % Min. Eps (is used for stabillity of the code as: dt<dr/2*max(velocity) for the stability)
    % (-1:Auto) (relative dielectric permitivity): Min. Epsilon to compute
    % Courant Stability Criteria:
    try
        disp(sprintf('%s %f','Min. eps. value (in Eps0):',addpar.Epsmin))
    catch
        addpar.Epsmin=1;
        disp(sprintf('%s %i %s','Min. eps. value (in Eps0):',addpar.Epsmin,'(default)'))
    end

    disp('-----------------------------------------------')
    % Source-Type (1:G, 2:DG, 3:R, 4:Ex):
    try
        switch addpar.srcWType
            case{1}
                disp(sprintf('%s %s','Source-wavelet-type:','Gauss'))
            case{2}
                disp(sprintf('%s %s','Source-wavelet-type:','Diff. Gauss'))
            case{3}
                disp(sprintf('%s %s','Source-wavelet-type:','Ricker'))
            case{4}
                disp(sprintf('%s %s','Source-wavelet-type:','External'))
        end
    catch
        addpar.srcWType=1;
        disp(sprintf('%s %s %s','Source-wavelet-type:','Gauss','(default)'))
    end

    % Inc/Exc Field at Source: Source implementation, either 1 (hard-source) or
    % 2 (soft-source, that accounts for the field at source-position):
    try
        switch addpar.srcImpl
            case{0}
                disp(sprintf('%s %s','Source implementation:','Hard-source'))
            case{1}
                disp(sprintf('%s %s','Source implementation:','Soft-source'))
        end
    catch
        addpar.srcImpl=1;
        disp(sprintf('%s %s %s','Source implementation:','Soft-source','(default)'))
    end

    % Source orientation (0:Ez, 1:Er/Ex):
    try
        switch addpar.srcOri
            case{0}
                disp(sprintf('%s %s','Source orientation:','Ez'))
            case{1}
                disp(sprintf('%s %s','Source implementation:','Er/Ex'))
        end
    catch
        addpar.srcOri=0;
        disp(sprintf('%s %s %s','Source orientation:','Ez','(default)'))
    end

    % Total Pulse Length in (s) if Source-type=1 or 2:
    % Center frequency if Source-type=3:
    try
        switch addpar.srcWType
            case{1,2}
                disp(sprintf('%s %e','Total pulse length (s):',addpar.Tg))
            case{3}
                disp(sprintf('%s %e','Center frequency (Hz):',addpar.Tg))
            otherwise
        end
    catch
        switch addpar.srcWType
            case{3}
                addpar.Tg=14*10^6;
                disp(sprintf('%s %e %s','Center frequency (Hz):',addpar.Tg,'(default)'))
            case{1,2}
                addpar.Tg=1.4*10^-8;
                disp(sprintf('%s %e %s','Total pulse length (s):',addpar.Tg,'(default)'))
            otherwise
        end
    end

    % Pulse Width only if Source-type=1 or 2:
    try
        switch addpar.srcWType
            case{1,2}
                disp(sprintf('%s %e','Pulse width (in s):',addpar.Tp))
            otherwise
        end
    catch
        switch addpar.srcWType
            case{1,2}
                addpar.Tp=1.6*10^-9;
                disp(sprintf('%s %e %s','Pulse width (in s):',addpar.Tp,'(default)'))
            otherwise
        end
    end

    disp('-----------------------------------------------')

    % Magnetic Permeability:
    try
        disp(sprintf('%s %f','Magnetic permeability (in Mu0):',addpar.mu))
    catch
        addpar.mu=1;
        disp(sprintf('%s %i %s','Magnetic permeability (in Mu0):',addpar.mu,'(default)'))
    end

    % Equivalent Magnetic Loss:
    try
        disp(sprintf('%s %f','Equivalent magnetic loss:',addpar.rho))
    catch
        addpar.rho=0;
        disp(sprintf('%s %i %s','Equivalent magnetic loss:',addpar.rho,'(default)'))
    end
    disp('-----------------------------------------------')

    % Output the size of the modes to the screen:
    disp(sprintf('%s %3.3f %s%i %s','Cell size:',addpar.dx,'m (',addpar.dx*100,'cm)'))
    disp(sprintf('%s %3.2f %s%i %s','Model width:',addpar.dx*addpar.Epsx,'m (',addpar.Epsx,'cells)'))
    disp(sprintf('%s %3.2f %s%i %s','Model depth:',addpar.dx*addpar.Epsz,'m (',addpar.Epsz,'cells)'))
    disp(sprintf('%s %3.2f %s','Simulation time:',addpar.t*10^9,'ns'))
    disp('-----------------------------------------------')
    try
        disp(sprintf('%s %3.3f %s%i %s','Inversion cell size:',addpar.dx_inv,'m (',addpar.dx_inv*100,'cm)'))
    catch
        addpar.dx_inv=0.1;
        disp(sprintf('%s %3.3f %s%i %s','Inversion cell size:',addpar.dx_inv,'m (',addpar.dx_inv*100,'cm)'))
    end
    
    try
        switch addpar.output
            case{1}
                disp('Output field: Ez')
            case{2}
                disp('Output field: Ex/Er')
            case{3}
                disp('Output field: Ez and Ex/Er')
            case{4}
                disp('Output field: Ez - file numbering defined by addpar.no')
            otherwise
                disp('Invalid choice of output')
                error=1;
                return
        end
    catch
        addpar.output=1;
        disp('Output field: Ez (default)')
    end
       
    try
        switch addpar.sim_mode
            case{1}
                disp('Simulaton mode: Pure forward')
            case{2}
                disp('Simulaton mode: Inversion: Forward-backward')
            case{3}
                disp('Simulaton mode: Inversion: Forward')
            otherwise
                disp('Invalid choice of output')
                error=1;
                return
        end
    catch
        addpar.sim_mode=1;
        disp('Simulation mode: Pure forward (default)')
    end
    
    % Indicat if this is the first iteration of multiple iteration:
    try
        if addpar.first==0
            disp('Using old parameter files')
        end
    catch
        addpar.first=1;
        disp('First iteration - parameter files are initialized')
    end
    
    disp('-----------------------------------------------')

    try
        disp(sprintf('%s %i','Number of cores applied:',addpar.cores))
    catch
        addpar.cores=1;
        disp(sprintf('%s %i %s','Number of cores applied:',addpar.cores,'(default)'))
    end

    % ========= Plot the models and receiver and transmitter setup ===========%
    try
        if addpar.plot==1
            figure_focus(1)
            Eps0    = 8.85418781762039080*1e-12; %in As/Vm = C^2/(N*m^2)
            subplot(121)
            imagesc(Eps./Eps0);
            cb1 = colorbar('horiz');
            set(get(cb1,'XLabel'),'String','Eps [Eps0]')
            hold on;
            plot(ant_pos(:,3)./addpar.dx,ant_pos(:,4)./addpar.dx,'rs');
            plot(ant_pos(:,3)./addpar.dx,ant_pos(:,4)./addpar.dx,'rx');
            plot(ant_pos(:,1)./addpar.dx,ant_pos(:,2)./addpar.dx,'go');
            plot(ant_pos(:,1)./addpar.dx,ant_pos(:,2)./addpar.dx,'gx');
            hold off;
            title('Dielectric permitivity (Eps)');
            xlabel('Width [cells]');
            ylabel(sprintf('Depth [cells] - (1 cell = %3.1fcm)',addpar.dx*100));
            axis equal;
            axis tight;
            subplot(122)
            imagesc(Sig*1000);
            cb2 = colorbar('horiz');
            set(get(cb2,'XLabel'),'String','Sig [mS/m]')
            hold on;
            plot(ant_pos(:,3)./addpar.dx,ant_pos(:,4)./addpar.dx,'rs');
            plot(ant_pos(:,3)./addpar.dx,ant_pos(:,4)./addpar.dx,'rx');
            plot(ant_pos(:,1)./addpar.dx,ant_pos(:,2)./addpar.dx,'go');
            plot(ant_pos(:,1)./addpar.dx,ant_pos(:,2)./addpar.dx,'gx');
            hold off;
            title('Electric conductivity (Sig)');
            xlabel('Width [cells]');
            ylabel(sprintf('Depth [cells] - (1 cell = %3.1fcm)',addpar.dx*100));
            axis equal;
            axis tight;
        else
            disp(sprintf('%s','Plot of the setup is deselected'))
            disp('-----------------------------------------------')
        end
    catch
        disp(sprintf('%s','Plot of the setup is deselected (default)'))
        disp('-----------------------------------------------')
    end

    % Check if the start-flag has been set to "on"
    try
        if addpar.start==1
            disp('Start of FDTD simulation: "on"')
        elseif addpar.start==0
            disp('Start of FDTD simulation: "off"')
        else
            disp('Warning: The start indicator value is incorrect')
            return
            error=1;
        end
    catch
        addpar.start=1;
        disp('Start of FDTD simulation: "on" (default)')
    end
    
    % Check if processing status is activated:
    try
        if addpar.status==1
            disp('Processing status screendump: "on"')
        elseif addpar.status==0
            disp('Processing status screendump: "off"')
        else
            disp('Warning: The screendump indicator value is incorrect')
            return
            error=1;
        end
    catch
        addpar.status=0;
        disp('Processing status screendump: "off" (default)')
    end

catch
    error=1;
    return
end
