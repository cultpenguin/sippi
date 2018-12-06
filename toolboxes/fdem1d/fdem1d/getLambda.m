function lambda = getLambda(r,j0len,j1len)
%% load Guptasarma1997 filters - make sure digFilt/ is in the path!
% B. Minsley, March 2010

%% JO
if (strcmp(j0len,'long'))
    load J0_120pt.mat;
    flen = 120;
else
    load J0_61pt.mat;
    flen = 61;
end
    
lambda.j0.a = a;
lambda.j0.s = s;
lambda.j0.w = w;
lambda.j0.flen = flen;

% make lambda for all coil spacings: Guptasarma1997, eqn2.
% lambda is size #coil-spacings (rows) by filter length (cols)
i = 1:flen;
l1 = 1./r;
l2 = 10.^(a+(i-1)*s);
lambda.j0.lam = l1 * l2;

%% J1
if (strcmp(j1len,'long'))
    load J1_140pt.mat;
    flen = 140;
else
    load J1_47pt.mat;
    flen = 47;
end
    
lambda.j1.a = a;
lambda.j1.s = s;
lambda.j1.w = w;
lambda.j1.flen = flen;
    

% make lambda for all coil spacings: Guptasarma1997, eqn2.
% lambda is size #coil-spacings (rows) by filter length (cols)
i = 1:flen;
l1 = 1./r;
l2 = 10.^(a+(i-1)*s);
lambda.j1.lam = l1 * l2;
