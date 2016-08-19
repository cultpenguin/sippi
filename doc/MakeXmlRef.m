i=0;
i=i+1;f{i}='.';
i=i+1;f{i}=['toolboxes',filesep,'gpr_fd'];
i=i+1;f{i}=['toolboxes',filesep,'frequency_matching'];
i=i+1;f{i}=['toolboxes',filesep,'traveltime'];
i=i+1;f{i}=['toolboxes',filesep,'covariance_inference'];
i=i+1;f{i}=['plotting'];
i=i+1;f{i}=['misc'];
%i=i+1;f{i}='sgems/*.m';
%i=i+1;f{i}='snesim/*.m';
%i=i+1;f{i}='misc/*.m';
%i=i+1;f{i}='fast/*.m';


for ff=1:length(f)
    
    FILES=dir(['..',filesep,f{ff},filesep,'*.m']);
    if ff==1
        name='sippi_functions.xml';
    
    else
        name=sprintf('sippi_functions_%s.xml',space2char(f{ff},'_',filesep));

    end
    fid=fopen(name,'w');
    
    for i=1:length(FILES)
      [p,name,ext]=fileparts(FILES(i).name);
      disp(name)
      h=help(name);
      
      %fprintf(fid,'<sect2 id=\"%s\"><title>%s</title>\n',name,name);
      fprintf(fid,'<sect2 xml:id=\"%s\"><title>%s</title>\n',name,name);
      fprintf(fid,'<para><programlisting><![CDATA[%s]]></programlisting></para>\n',h);
      fprintf(fid,'</sect2>\n\n');
      
    end

    try
        fclose(fid);
    end
    
end
  
