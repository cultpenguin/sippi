function [prd, varargout] = fdem1dfwd(S,M,el,flg)
% given a system, model, and coordinates, predict data
% optionally output sensitivities if flg = 1;
% B. Minsley, June 2010
%   *** changed C.el to el, May 2012
%
% From 
% Elwaseif, M., Robinson, J., Day-Lewis, F. D., Ntarlagiannis, D., Slater, L. D., Lane, J. W., ... & Schultz, G. (2017). A matlab-based frequency-domain electromagnetic inversion code (FEMIC) with graphical user interface. Computers & Geosciences, 99, 61-71.
% Minsley, B. J. (2011). A trans-dimensional Bayesian Markov chain Monte Carlo algorithm for model assessment using frequency-domain electromagnetic data. Geophysical Journal International, 187(1), 252-272.


%% set lambda values as global variable (once per system configuration)
persistent lambda
if isempty(lambda)
    lambda = getLambda(S.r,'long','long');
end

%% switch for coil orientations
zz = find(strcmp(S.tor,'z') & strcmp(S.ror,'z')); %HCP
yy = find(strcmp(S.tor,'y') & strcmp(S.ror,'y')); %VCX
xx = find(strcmp(S.tor,'x') & strcmp(S.ror,'x')); %VCP

yz = find(strcmp(S.tor,'z') & strcmp(S.ror,'y')); %TZ-RY
xz = find(strcmp(S.tor,'z') & strcmp(S.ror,'x')); %TZ-RX

zx = find(strcmp(S.tor,'x') & strcmp(S.ror,'z')); %TX-RZ
yx = find(strcmp(S.tor,'x') & strcmp(S.ror,'y')); %TX-RY

zy = find(strcmp(S.tor,'y') & strcmp(S.ror,'z')); %TY-RZ
xy = find(strcmp(S.tor,'y') & strcmp(S.ror,'x')); %TY-RX

if flg == 0  %% no sensitivity output
    %% calculate reflection coefficient
    %tic
    if ~isempty(zz) || ~isempty(xx)|| ~isempty(yy) || ~isempty(xy) || ~isempty(yx)
        [rTE.j0,u0.j0] = calcrTEsens(S,M,lambda.j0.lam,flg);
    end
    if ~isempty(zx) || ~isempty(zy) || ~isempty(xz) || ~isempty(yz) || ~isempty(xx) ...
            || ~isempty(yy) || ~isempty(xy) || ~isempty(yx)
        [rTE.j1,u0.j1] = calcrTEsens(S,M,lambda.j1.lam,flg);
    end
    %a=toc
    
    %% calculate H & H0 for all orientations
    % z-oriented transmitter
    %tic
    if ~isempty(zz) %HCP
        [H(zz),H0(zz)] = calcHzz(zz,S,el,rTE.j0,u0.j0,lambda.j0);
    end
    if ~isempty(xx) %VCX
        [H(xx),H0(xx)] = calcHxx(xx,S,el,rTE,u0,lambda);
    end
    %b=toc
    %[a/(a+b) b/(a+b)]
    % if ~isempty(yz)
    %     H0(yz) = calcHyz(yz,S,1e3*C.z(1),rTE.j1,u0.j1,lambda.j1);
    %     H(yz) = calcHyz(yz,S,C.z(1),rTE.j1,u0.j1,lambda.j1);
    % end
    % if ~isempty(xz)
    %     H0(xz) = calcHxz(xz,S,1e6*C.z(1),rTE.j1,u0.j1,lambda.j1);
    %     H(xz) = calcHxz(xz,S,C.z(1),rTE.j1,u0.j1,lambda.j1);
    % end
    %
    % % x-oriented transmitter
    % if ~isempty(zx)
    %     H0(zx) = calcHzx(zx,S,1e6*C.z(1),rTE.j1,u0.j1,lambda.j1);
    %     H(zx) = calcHzx(zx,S,C.z(1),rTE.j1,u0.j1,lambda.j1);
    % end

    % scaling factor for moments- need to fix this...
    scl = S.tmom .* S.rmom;
    
    prd(:,1) = 1e6*(H - H0)./H0;
    prd = prd.*scl;
    
    varargout = {};

elseif flg == 1 %% output sensitivities
        %% calculate reflection coefficient
    if ~isempty(zz) || ~isempty(xx)|| ~isempty(yy) || ~isempty(xy) || ~isempty(yx)
        [rTE.j0,u0.j0,drTEdlogcon.j0] = calcrTEsens(S,M,lambda.j0.lam,flg);
    end
    if ~isempty(zx) || ~isempty(zy) || ~isempty(xz) || ~isempty(yz) || ~isempty(xx) ...
            || ~isempty(yy) || ~isempty(xy) || ~isempty(yx)
        [rTE.j1,u0.j1,drTEdlogcon.j1] = calcrTEsens(S,M,lambda.j1.lam,flg);
    end
    
    %% calculate H & H0 for all orientations plus sensitivities
    
       
    for i = 1:length(M.con) % loop over layers
        
        if ~isempty(zz) %HCP
            if i==1,[H(zz),H0(zz)] = calcHzz(zz,S,el,rTE.j0,u0.j0,lambda.j0);end
            [dH(zz,i),dH0(zz,i)] = calcHzz(zz,S,el,drTEdlogcon.j0(:,:,i),u0.j0,lambda.j0);
        end
        if ~isempty(xx) %VCX
            if i==1,[H(xx),H0(xx)] = calcHxx(xx,S,el,rTE,u0,lambda);end
            drTE.j0 = drTEdlogcon.j0(:,:,i);
            drTE.j1 = drTEdlogcon.j1(:,:,i);
            [dH(xx,i),dH0(xx,i)] = calcHxx(xx,S,el,drTE,u0,lambda);
        end
        % if ~isempty(yz)
        %     H0(yz) = calcHyz(yz,S,1e3*C.z(1),rTE.j1,u0.j1,lambda.j1);
        %     H(yz) = calcHyz(yz,S,C.z(1),rTE.j1,u0.j1,lambda.j1);
        % end
        % if ~isempty(xz)
        %     H0(xz) = calcHxz(xz,S,1e6*C.z(1),rTE.j1,u0.j1,lambda.j1);
        %     H(xz) = calcHxz(xz,S,C.z(1),rTE.j1,u0.j1,lambda.j1);
        % end
        %
        % % x-oriented transmitter
        % if ~isempty(zx)
        %     H0(zx) = calcHzx(zx,S,1e6*C.z(1),rTE.j1,u0.j1,lambda.j1);
        %     H(zx) = calcHzx(zx,S,C.z(1),rTE.j1,u0.j1,lambda.j1);
        % end

    
    end
    
    % scaling factor for moments- need to fix this...
    scl = S.tmom .* S.rmom;
    
    prd(:,1) = 1e6*(H - H0)./H0;
    prd = prd.*scl;
    J = 1e6*(dH - dH0)./dH0;
    
    % Changed code here: ,Akbar
    %J = J.*repmat(scl,1,size(J,2)); 
    J = J.* scl(:,ones(1,size(J,2)));
    
    varargout = {J};

end