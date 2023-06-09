function write_workbench_xyz(fname,D,HEAD,TXT, nanval, NL,sep_type);

if nargin<2
    help(mfilename);
    return;
end

[nd,nc]=size(D);

if nargin<3;
    for ic=1:nc;
        HEAD{ic}=sprintf('DATA_%d',ic);
    end
end

if nargin<4;
    TXT{1}='Matlab Export in WorkBench XYZ format';
end

if nargin<5;
    nanval = -9999;
end

if nargin<6;
    NL = [];
end

if nargin<7;
    sep_type = 0; % ''
    %sep_type = 1;% ';'
    %sep_type = 2;% ',' 
end


if ~isempty(NL);
    i=length(TXT);
    TXT{i+1}=sprintf('NUMBER OF LAYERS');
    TXT{i+2}=sprintf('%d',NL);
end
if ~isempty(nanval);
    D(isnan(D))=nanval;
end


%nchar = 14; % number of chars per columns
disp(sprintf('%s: Writing XYZ file to ''%s''',mfilename,fname))
%% START IO
fid = fopen(fname,'w');

% write top lines
for i=1:length(TXT);
    if length(TXT{i})==0,TXT{i}=' ';end 
    if ~strcmp(TXT{i}(1),'/');
        TXT{i}=['/',TXT{i}];
    end
    %fprintf(fid,'%s\r\n',TXT{i});
    fprintf(fid,'%s\r\n',TXT{i});
    
     
end

% write headers
fprintf(fid,'/%13s  ',HEAD{1});
for ic=2:nc;
    fprintf(fid,'%14s',HEAD{ic});
end
fprintf(fid,'\r\n');

for id=1:nd;
    if sep_type==0;
        fprintf(fid,' %13f',D(id,1:(end-1)));
    elseif sep_type==1;
        fprintf(fid,' %12f;',D(id,1:(end-1)));
    else
        fprintf(fid,' %12f,',D(id,1:(end-1)));
    end
    fprintf(fid,' %12f',D(id,end));

    fprintf(fid,'\r\n');
end

fclose(fid);


