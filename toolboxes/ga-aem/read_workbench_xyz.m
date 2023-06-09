function [D,HEAD,TXT,GATE_TIME]=read_workbench_xyz(filename,l_head,nanval)

if nargin < 1, filename='IGIS_export_MOD_dat.xyz';end

disp(sprintf('%s: reading %s',mfilename,filename))
fid = fopen(filename,'r');

if nargin < 2,
    % try get the number of lines in the geader
    header = textscan(fid, '%s', 100, 'Delimiter', '\n');
    for il=1:length(header{1});
        if strcmp(header{1}{il}(1),'/')==0
            l_head = il-1;
            break
        end
    end
    fseek(fid, 0, 'bof');
end
if nargin < 3, nanval = 99999;end
%clear all
%filename='IGIS_export_MOD_dat.xyz';
%l_head = 7;



for i=1:(l_head-1);
    TXT{i} = fgetl(fid);
end


headline = fgetl(fid);
HEAD = split(headline);
HEAD = HEAD(2:end);
nc = length(HEAD);

f_position_data = ftell(fid);

allReals=1;
%% First try reading all data as floats.
if allReals == 1
    try
        D = fscanf(fid,'%f');
        nd=length(D);
        nr = nd/nc;
        D=reshape(D,[nc,nr])';
    catch
        disp(sprintf('%s: Could not read %s as all floats. trying as chars',mfilename,filename))
        allReals=0;
    end
end

%% MIXED TYPE
if allReals==0;
    fseek(fid, f_position_data, 'bof');
    Dchar = textscan(fid, '%s');
    nd=length(Dchar{1});
    nrows = nd/nc;

    D = zeros(nrows,nc);
    for ir=1:nrows
        if mod(ir,100)==0,progress_txt(ir,nrows,filename);end
        for ic=1:nc
            i = (ir-1)*nc+(ic-1)+1;
            try;D(ir,ic) = str2num(Dchar{1}{i});end
        end
    end
end


if ~isempty(nanval);
    D(D==nanval)=NaN;
end

%% CHECK FOR GATE TIMES
igates=[];
igatetime=[];
for i=1:length(TXT);
    if strcmp(upper(TXT{i}),'/NUMBER OF GATES');
        igates = i+1;
    end
    if strfind(upper(TXT{i}),'/GATE TIMES');
        igatetime = i+1;
    end
end

if ~isempty(igates)
    NGATES = str2num(TXT{igates}(2:end));
else
    NGATES=[];
end
if ~isempty(igatetime)
    GATE_TIME = str2num(TXT{igatetime}(2:end));
else
    GATE_TIME=[];
end

disp(sprintf('%s: read %s, nc=%d, nr=%d',mfilename,filename,nc,size(D,1)));

fclose(fid);