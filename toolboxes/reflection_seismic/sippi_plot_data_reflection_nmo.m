% sippi_plot_data_reflection_nmo: Plot data response for sippi_forward_reflection_nmo
% Call:
%    sippi_plot_data_reflection_nmo(d,data,i_data,prior);
%
%
% See also sippi_plot_data, sippi_forward_reflection_nmo
%
function sippi_plot_data(d,data,id_arr,prior);


if nargin<3;
  id_arr=1:length(d);
end

for id=id_arr;
  
  if ~ishold
    figure_focus(20+id);
  end
  hold off
  if nargin>3
    % wiggle plot
    nt=length(prior{1}.y);
    t=prior{1}.y;
    n=length(d{id});
    na=n/nt;
    
    
    if ~isempty(data)
      %hold on
      avo_gather_obs=reshape(data{id}.d_obs,nt,na);
      dmax=wiggle(1:na,t,avo_gather_obs,'VA');
      %hold off
    end
    
    hold on
    avo_gather=reshape(d{id},nt,na);
    wiggle(1:na,t,avo_gather,'wiggle',dmax);
    hold off
  else
    % simple plot
    plot(d{id},'k-')
    if nargin>1
      hold on
      plot(data{id}.d_obs,'r-')
    end
    
  end
  
end
hold off

