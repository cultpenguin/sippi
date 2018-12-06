function D = readData(dataFile,nf)
% expected input format as .csv (including header line)
% EL is sensor elevation (positive downward)
% nf is number of frequencies
% LINE FID  X  Y  Z  EL  I1  Q1 ... In  Qn I1err Q1err ... Inerr Qnerr

dat = importdata(dataFile,',',1);
nrec = size(dat.data,1);

D = struct([]);
for i = 1:nrec
    D(i).line = dat.data(i,1);
    D(i).fid = dat.data(i,2);
    D(i).x = dat.data(i,3);
    D(i).y = dat.data(i,4);
    D(i).z = dat.data(i,5);
    D(i).el = dat.data(i,6);
    D(i).obs = dat.data(i,7:7+2*nf-1)';
    D(i).err = dat.data(i,7+2*nf:end)';
end
    