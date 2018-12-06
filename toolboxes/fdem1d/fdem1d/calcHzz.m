% Modified by Akbar
function [Hzz,Hzz0] = calcHzz(ix,S,z,rTE,u0,lambda);
%TZ-RZ
% B. Minsley, March 2010

%f(r) = int(K(lam)*Ji(lam*r)dlam
%r*f(r) = sum(K(lam)*W)
 
% decompose
w = lambda.w(:);
flen = lambda.flen;
lam = lambda.lam(ix,:);

% z is positive downwards 
h = -(z + S.tzoff(ix)); % transmitter height
rz = z + S.rzoff(ix);   % receiver elevation
% Function that changed can be seen as text, Akbar

%H = repmat(h,1,flen);
H = h(:,ones(1,flen));

%Z = repmat(rz,1,flen);
Z = rz(:,ones(1,flen));


% truncate
u0 = u0(ix,:);
rTE = rTE(ix,:);
tmom = S.tmom(ix);


% K(lam) - eqn. 4.46
K = ( exp(-u0.*(Z+H)) + rTE .* exp(u0.*(Z-H)) ) .* (lam.^3) ./ u0;
% K = ( exp(-u0.*(rzoff)) + rTE.* exp(u0.*(2*tz+rzoff)) ) .* lam(.^3) ./ u0;

Hzz = (tmom/(4*pi)) .* (K*w) ./ S.r(ix);

% free-space
K0 = ( exp(-u0.*(Z+H)) ) .* (lam.^3) ./ u0;

Hzz0 = (tmom/(4*pi)) .* (K0*w) ./ S.r(ix);