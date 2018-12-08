function [prd,M] = calcFwd(S,model,htx,ds,i)
% depth grid to calculate forward response- careful it is not too coarse 
% otherwise forward responses will be inaccurate!
%   for the current example, this downsampling leads to <3% error, ~15ppm
%   max, which is acceptable
% B. Minsley Nov 2013

if ds
    zfwd = [(0:9)'; logspace(log10(10),log10(max(model.Z(:))),15)'];
    zbnd = [model.Z(:,1) [model.Z(2:end,1);10*model.Z(end,1)]];
else
    zfwd = model.Z(:,1);
end


% downsample model (if flagged) and calculate forward response.
M.k = length(zfwd);
M.z = zfwd;
M.thk = [diff(M.z);0];
M.chie = zeros(M.k,1);
M.chim = zeros(M.k,1);
if ds
    for j = 1:length(zfwd)
        if j < length(zfwd)
            thkin = max(0,min(zbnd(:,2),zfwd(j+1)) - max(zbnd(:,1),zfwd(j)));
        else
            thkin = max(0,min(zbnd(:,2),10*zfwd(j)) - max(zbnd(:,1),zfwd(j)));
        end
        M.con(j,1) = thkin' * model.CON(:,i) / sum(thkin);
    end
else
    M.con = model.CON(:,i);
end

%calculate and store forward response, separating real and imaginary
%components
prdcmplx = fdem1dfwd(S,M,htx,0);
prd(1:2:12,1)=real(prdcmplx);
prd(2:2:12,1)=imag(prdcmplx);