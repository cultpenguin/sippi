% Modified by Akbar
function [Hxx,H0xx] = calcHxx(ix,S,z,rTE,u0,lambda);
%%VCX
% B. Minsley, March 2010

%f(r) = int(K(lam)*Ji(lam*r)dlam
%r*f(r) = sum(K(lam)*W)

% decompose
w.j0 = lambda.j0.w(:);
flen.j0 = lambda.j0.flen;
lam.j0 = lambda.j0.lam(ix,:);

w.j1 = lambda.j1.w(:);
flen.j1 = lambda.j1.flen;
lam.j1 = lambda.j1.lam(ix,:);

% z is positive downwards 
h = -(z + S.tzoff(ix)); % transmitter height
rz = z + S.rzoff(ix);   % receiver elevation
% Function that changed can be seen as text, Akbar

%H.j0 = repmat(h,1,flen.j0);
H.j0 = h(:,ones(1,flen.j0));

%Z.j0 = repmat(rz,1,flen.j0);
Z.j0 = rz(:,ones(1,flen.j0));

%H.j1 = repmat(h,1,flen.j1);
H.j1 = h(:,ones(1,flen.j1));

%Z.j1 = repmat(rz,1,flen.j1);
Z.j1 = rz(:,ones(1,flen.j1));

% truncate
u0.j0 = u0.j0(ix,:);
rTE.j0 = rTE.j0(ix,:);
tmom = S.tmom(ix);

u0.j1 = u0.j1(ix,:);
rTE.j1 = rTE.j1(ix,:);


% K(lam) - eqn. 4.46
K.j1 = ( exp(-lam.j1.*(Z.j1+H.j1)) - rTE.j1 .* exp(lam.j1.*(Z.j1-H.j1)) ) .* (lam.j1);
K.j0 = ( exp(-lam.j0.*(Z.j0+H.j0)) - rTE.j0 .* exp(lam.j0.*(Z.j0-H.j0)) ) .* (lam.j0).^2;

Hxxj1 = -(tmom/(4*pi)) .* ( (1./S.r(ix)) - (2*S.rx(ix).^2./S.r(ix).^3) ) .* (K.j1*w.j1) ./ S.r(ix);
Hxxj0 = -(tmom/(4*pi)) .* ( S.rx(ix)./S.r(ix) ).^2 .* (K.j0*w.j0) ./ S.r(ix);

Hxx = Hxxj1 + Hxxj0;

% free-space
K0.j1 = ( exp(-lam.j1.*(Z.j1+H.j1)) ) .* (lam.j1) ;
K0.j0 = ( exp(-lam.j0.*(Z.j0+H.j0)) ) .* (lam.j0).^2 ;

H0xxj1 = -(tmom/(4*pi)) .* ( (1./S.r(ix)) - (2*S.rx(ix).^2./S.r(ix).^3) ) .* (K0.j1*w.j1) ./ S.r(ix);
H0xxj0 = -(tmom/(4*pi)) .* ( S.rx(ix)./S.r(ix) ).^2 .* (K0.j0*w.j0) ./ S.r(ix);

H0xx = H0xxj1 + H0xxj0;

