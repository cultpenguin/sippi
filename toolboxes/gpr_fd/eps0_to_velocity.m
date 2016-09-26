

function eps0_to_velocity(input)

c0=2.99792458*10^8;
velo=sqrt(c0.^2/input); % input=EPS0

disp([])
if velo*10^-9<0.3
disp(sprintf('%2.3f eps/eps0 ~ %2.3f m/ns',input,velo*10^-9))
end
disp([])
EPS0=c0^2/(input*10^9)^2; %input=velo
if EPS0>1
disp(sprintf('%2.3f m/ns ~ %2.3f eps/eps0',input,EPS0))
end