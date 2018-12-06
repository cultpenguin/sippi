function [] = oneDmod(z,v,col)

% input depths & values should have same length
% each value is associated with a top depth given by the corresponding z

z=z(:);
v=v(:);

plot(kron(v,[1;1]),[0;kron(z,[1;1]);1.5*max(abs(z))*sign(max(z))],[col '-'])

