% Modified by Akbar
% Function that changed were too numerous to keep in
function [varargout] = calcrTEsens(S,M,lam,flg)
%% calculate reflection coefficient according to Ward and Hohmann, EM theory for geophysical applications
% B. Minsley, June 2010

% constants
eps0 = 1/(35950207149.4727056*pi);%8.8541878176e-12;
mu0 = 4*pi*1e-7;

flen = size(lam,2);

% angular frequency
omega = 2*pi*S.freq;
nf = length(omega);

%% add air layer onto top of model
con = [0 M.con'];
chie = [0 M.chie'];
chim = [0 M.chim'];
nlyr = length(con);
h = [NaN;M.thk];

%% set impedivity, admitivity, wavenumber (zn,yn,un)
% un = sqrt(kx^2 + ky^2 - kn^2) = sqrt(lam^2 - kn^2)
% kn = sqrt(-zn*yn) = sqrt(omega^2*mu_n*eps_n - j*omega*mu_n*sigma_n)
% zn = j*omega*mu_n, yn = j*omega*eps_n + sigma_n
yn = 1j*omega*eps0*(1+chie) + con(ones(nf,1),:);
yn = reshape(yn,nf,1,nlyr);
zn = 1j*omega*mu0*(1+chim);
zn = reshape(zn,nf,1,nlyr);

% eqn 4.29
lam_sq = lam.^2;
un = sqrt( lam_sq(:,:,ones(1,1,nlyr)) + ...
    zn(:,ones(1,flen,1),:).*yn(:,ones(1, flen, 1),:) );

%% set layer admittances, Yn
% eqn 4.27
Yn = un ./ zn(:,ones(1, flen, 1),:);

if flg == 0 
    %% no sensitivity output
    % set half-space admittance derivative wrt log(sigma)
    % eqn 4.26 (Y is Yhat)
    Y(:,:,nlyr) = Yn(:,:,nlyr);
    Y(:,:,1) = NaN; % air layer
    
    % in case of a half-space
    if nlyr == 2
        rTE = (Yn(:,:,1) - Y(:,:,2)) ./ (Yn(:,:,1) + Y(:,:,2));
        u0 = un(:,:,1);
        varargout = {rTE,u0};
        return
    end
    
    %% iterate to find surface admittance & sensitivity
    for i = nlyr-1:-1:2
        tanuh = tanh(un(:,:,i)*h(i)); % store common operation
        
        % Y = Yhat for each layer
        num = Y(:,:,i+1) + Yn(:,:,i) .* tanuh;
        den = Yn(:,:,i) + Y(:,:,i+1) .* tanuh;
        Y(:,:,i) = Yn(:,:,i) .* num ./ den;
    end
    
    % outputs
    rTE = (Yn(:,:,1) - Y(:,:,2)) ./ (Yn(:,:,1) + Y(:,:,2));
    u0 = un(:,:,1);
    
    varargout = {rTE,u0};
    
else
    %% initialize arrays for admittance & sensitivities
    % set half-space admittance derivative wrt log(sigma)
    % eqn 4.26 (Y is Yhat)
    Y(:,:,nlyr) = Yn(:,:,nlyr);
    Y(:,:,1) = NaN; % air layer
    dYdlogcon(:,:,nlyr) = con(end)./(2*un(:,:,end));
    dYdlogcon(:,:,1) = NaN; % air layer
    dYdY = ones([size(lam) nlyr]);
    
    % in case of a half-space
    if nlyr == 2
        rTE = (Yn(:,:,1) - Y(:,:,2)) ./ (Yn(:,:,1) + Y(:,:,2));
        u0 = un(:,:,1);
        drTEdlogcon = -2*Yn(:,:,1).*dYdlogcon(:,:,nlyr)./(Yn(:,:,1) + Y(:,:,2)).^2;
        varargout = {rTE,u0,drTEdlogcon};
        return
    end
    
    %% iterate to find surface admittance & sensitivity
    for i = nlyr-1:-1:2
        tanuh = tanh(un(:,:,i)*h(i)); % store common operation
        
        % Y = Yhat for each layer
        num = Y(:,:,i+1) + Yn(:,:,i) .* tanuh;
        den = Yn(:,:,i) + Y(:,:,i+1) .* tanuh;
        Y(:,:,i) = Yn(:,:,i) .* num ./ den;
        
        % dY(i) / d(ln con(i))
        dYdlogcon(:,:,i) = (con(i)./(2*un(:,:,i).*den.^2)) .* (...
            2*Yn(:,:,i) .* Y(:,:,i+1) .* tanuh.^2 + ...
           (1j*omega(:,ones(1,flen))*mu0*(1+chim(i))*h(i)) .* (Yn(:,:,i).*Y(:,:,i+1).^2 - Yn(:,:,i).^3) .* tanuh.^2 + ...
           (Y(:,:,i+1).^2 - Yn(:,:,i).^2) .* tanuh + 2*Yn(:,:,i).^2 + ...
           (1j*omega(:,ones(1,flen))*mu0*(1+chim(i))*h(i)) .* (Yn(:,:,i).^3 - Yn(:,:,i).*Y(:,:,i+1).^2) );
       
        % dY(i) / dY(i+1)
        numYY = Yn(:,:,i).^2 .* (1 - tanuh.^2);
        denYY = ( Yn(:,:,i) + Y(:,:,i+1) .* tanuh ).^2;
        dYdY(:,:,i+1) = numYY ./ denYY;
        
    end
    dYdlogcon(:,:,1) = []; %air layer
    dYdY(:,:,1) = []; %air layer
    
    % outputs
    rTE = (Yn(:,:,1) - Y(:,:,2)) ./ (Yn(:,:,1) + Y(:,:,2));
    u0 = un(:,:,1);
    
    num = -2 * repmat(Yn(:,:,1),[1 1 nlyr-1]);
    den = repmat( (Yn(:,:,1) + Y(:,:,2)) .^2 , [1 1 nlyr-1] );
    
    drTEdlogcon = num .* (cumprod(dYdY,3) .* dYdlogcon) ./ den;
    
    varargout = {rTE,u0,drTEdlogcon};
end

