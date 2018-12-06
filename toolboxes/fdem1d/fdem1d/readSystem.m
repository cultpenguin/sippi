function S = readSystem(systemFile)

fid = fopen(systemFile);

sys = textscan(fid,'%f %s %f %f %f %f %s %f %f %f %f');
% sys file has the following columns:
% freq tor tmom tx ty th ror rmom rx ry rh 
% everything is relative to the system, i.e. tx/ty/th are relative to 
%  the survey x/y/z coordinates, and rx/ry/rh are offset from tx/ty/th
%  z is positive downwards 
S.freq = sys{1};
S.tor = sys{2};
S.tmom = sys{3};
S.tx = sys{4};
S.ty = sys{5};
S.tzoff = sys{6};
S.ror = sys{7};
S.rmom = sys{8};
S.rx = sys{9};
S.ry = sys{10};
S.rzoff = sys{11};

S.nf = length(S.freq);
S.r = sqrt( (S.tx - S.rx).^2 + (S.ty - S.ry).^2 );

fclose(fid);
