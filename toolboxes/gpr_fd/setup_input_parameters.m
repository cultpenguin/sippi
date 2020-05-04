
function [addpar error]=setup_input_parameters(ant_pos,Sig,Eps,addpar)
error=0;

try
    %=============== Parameters used in the forward modelling ================%


    % Name of forward controler #1:
    try
        isempty(addpar.forwardexe1);
    catch
        addpar.forwardexe1='FDTD_forward1.exe';
    end

    % Name of forward controler #2:
    try
        isempty(addpar.forwardexe2);
    catch
        addpar.forwardexe2='FDTD_forward2.exe';
    end

    % Name of forward controler #3:
    try
        isempty(addpar.forwardexe3);
    catch
        addpar.forwardexe3='FDTD_forward3.exe';
    end

    % Name of forward controler #4:
    try
        isempty(addpar.forwardexe4);
    catch
        addpar.forwardexe4='FDTD_forward4.exe';
    end
    
    % Name of forward controler #5:
    try
        isempty(addpar.forwardexe5);
    catch
        addpar.forwardexe5='FDTD_forward5.exe';
    end

    % Name of forward controler #6:
    try
        isempty(addpar.forwardexe6);
    catch
        addpar.forwardexe6='FDTD_forward6.exe';
    end

    % Name of forward controler #7:
    try
        isempty(addpar.forwardexe7);
    catch
        addpar.forwardexe7='FDTD_forward7.exe';
    end

    % Name of forward controler #8:
    try
        isempty(addpar.forwardexe8);
    catch
        addpar.forwardexe8='FDTD_forward8.exe';
    end
    
    % Name of forward controler #9:
    try
        isempty(addpar.forwardexe9);
    catch
        addpar.forwardexe9='FDTD_forward9.exe';
    end

    % Name of forward controler #10:
    try
        isempty(addpar.forwardexe10);
    catch
        addpar.forwardexe10='FDTD_forward10.exe';
    end

    % Name of forward controler #11:
    try
        isempty(addpar.forwardexe11);
    catch
        addpar.forwardexe11='FDTD_forward11.exe';
    end

    % Name of forward controler #12:
    try
        isempty(addpar.forwardexe12);
    catch
        addpar.forwardexe12='FDTD_forward12.exe';
    end
    
    % Name of forward controler #13:
    try
        isempty(addpar.forwardexe13);
    catch
        addpar.forwardexe13='FDTD_forward13.exe';
    end

    % Name of forward controler #14:
    try
        isempty(addpar.forwardexe14);
    catch
        addpar.forwardexe14='FDTD_forward14.exe';
    end

    % Name of forward controler #15:
    try
        isempty(addpar.forwardexe15);
    catch
        addpar.forwardexe15='FDTD_forward15.exe';
    end

    % Name of forward controler #16:
    try
        isempty(addpar.forwardexe16);
    catch
        addpar.forwardexe16='FDTD_forward16.exe';
    end

    % Number of boundary cells:
    try
        isempty(addpar.bx);
    catch
        addpar.bx=40;
    end

    % GPML-Tuning: sm (-1:auto): Max conductivity:
    try
        isempty(addpar.sm);
    catch
        addpar.sm=-1;
    end

    % GPML-Tuning: sigm (-1:auto): Max coord. strecking:
    try
        isempty(addpar.sigm);
    catch
        addpar.sigm=-1;
    end

    % GPML-Tuning expo: Exponent, describing how fast the strecting and
    % conductivity happens in the boundaries:
    try
        isempty(addpar.n);
    catch
        addpar.n=2;
    end

    % C-Stability-Scaling (do not change): Stability modifier, controls
    % dispersion criteria (Courant):
    try
        isempty(addpar.sc);
    catch
        addpar.sc=0.9;
    end

    % Min. Eps (is used for stabillity of the code as: dt<dr/2*max(velocity) for the stability)
    % (-1:Auto) (relative dielectric permitivity): Min. Epsilon to compute
    % Courant Stability Criteria:
    try
        isempty(addpar.Epsmin);
    catch
        addpar.Epsmin=1;
    end

    % Source-Type (1:G, 2:DG, 3:R, 4:Ex):
    try
        isempty(addpar.srcWType);
    catch
        addpar.srcWType=1;
    end

    % Inc/Exc Field at Source: Source implementation, either 1 (hard-source) or
    % 2 (soft-source, that accounts for the field at source-position):
    try
        isempty(addpar.srcImpl);
    catch
        addpar.srcImpl=1; % Soft source
    end

    % Source orientation (0:Ez, 1:Er/Ex):
    try
        isempty(addpar.srcOri);
    catch
        addpar.srcOri=0;
    end

    % Total Pulse Length in (s) if Source-type=1 or 2:
    % Center frequency if Source-type=3:
    try
        switch addpar.srcWType
            case{1,2}
                isempty(addpar.Tg);

            case{3}
                isempty(addpar.Tg);
        end

    catch
        switch addpar.srcWType
            case{3}
                addpar.Tg=14*10^6;
            case{1,2}
                addpar.Tg=1.4*10^-8;
        end
    end

    % Pulse Width only if Source-type=1 or 2:
    try
        switch addpar.srcWType
            case{1,2}
                isempty(addpar.Tp);
        end
    catch
        switch addpar.srcWType
            case{1,2}
                addpar.Tp=1.6*10^-9;
        end
    end

    % Magnetic Permeability:
    try
        isempty(addpar.mu);
    catch
        addpar.mu=1;
    end

    % Equivalent Magnetic Loss:
    try
        isempty(addpar.rho);
    catch
        addpar.rho=0;
    end

    % Output the size of the modes to the screen:
    try
        isempty(addpar.dx_inv);
    catch
        addpar.dx_inv=0.1;
    end

    try
        isempty(addpar.output);

    catch
        addpar.output=1;
    end

    try
        isempty(addpar.sim_mode);
    catch
        addpar.sim_mode=1;
    end

    try
        isempty(addpar.cores);
    catch
        addpar.cores=1;
    end
    
    try
        isempty(addpar.first);
    catch
        addpar.first=1;
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
            title('Dielectirc permitivity (Eps)');
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

        end
    end

    % Check if the start-flag has been set to "on"
    try
        isempty(addpar.start);
    catch
        addpar.start=1;
    end

    % Check if processing status is activated:
    try
        isempty(addpar.status==1);
    catch
        addpar.status=0;
    end

catch
    error=1;
    return
end
