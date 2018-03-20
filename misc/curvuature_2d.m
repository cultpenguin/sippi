function curvature=curvuature_2d(x,y,dx,dy);
if nargin<3, dx=ones(size(x));end
if nargin<y, dx=ones(size(y));end

dx  = gradient(x);
ddx = gradient(dx);
dy  = gradient(y);
ddy = gradient(dy);
num   = dx .* ddy - ddx .* dy;
denom = dx .* dx + dy .* dy;
denom = sqrt(denom);
denom = denom .* denom .* denom;
curvature = num ./ denom;
curvature(denom < 0) = NaN;