% sippi_plot_data_reflection_nmo: Plot data response for sippi_forward_reflection_nmo
% Call:
%    sippi_plot_data_reflection_nmo(d,data,i_data,prior);
%
%
% See also sippi_plot_data, sippi_forward_reflection_nmo
%
function sippi_plot_data_reflection_nmo(d,data,id_arr,options);


if nargin<3;
  id_arr=1:length(d);
end
if nargin<4,options.null='';end

if isfield(options,'na'); 
    na=options.na;
else; 
    na=1;
end

for id=id_arr;
  
  if ~ishold
    figure_focus(20+id);
  end
  
  nt=length(d{id})/na;
  if isfield(options,'t'); 
      t=options.t;
      nt=length(t);
      na=length(d{id})/nt;
  else, 
      t=1:nt;
  end

  avo_gather=reshape(d{id},nt,na);
    
  
  
  hold off
  if nargin>3
    if ~isempty(data)
      %hold on
      avo_gather_obs=reshape(data{id}.d_obs,nt,na);
      dmax=wiggle(1:na,t,avo_gather_obs,'VA');
      %hold off
      hold on
      wiggle(1:na,t,avo_gather,'wiggle',dmax);
      hold off
    else
        wiggle(1:na,t,avo_gather,'wiggle');
    end
    
  else
       % simple plot
       wiggle(1:na,t,avo_gather,'VA');    
  end
  
end
hold off

