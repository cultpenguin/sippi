function id=xyz_get_header_id(H,str)
id=[];
for i=1:length(H)
    if strcmp(H{i},str)
        id=i;
    end
end