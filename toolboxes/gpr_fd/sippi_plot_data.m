% sippi_plot_data: plot GPR data
%                  overwrites sippi_plot_data plot data in SIPPI
%
% Call.
%    sippi_plot_data(d,data);
%
function sippi_plot_data(d,data);

figure_focus(20+1);clf;set_paper('portrait')    
for id=1:length(d);
    
  D1(:,id)=d{id};
  
  if nargin==2;
    D2(:,id)=data{id}.d_obs;
  end
 
end
[nt,nx]=size(D1);


wiggle(1:1:nx,1:1:nt,D1,'VA');
if nargin>1
  subplot(3,1,1);
  wiggle(1:1:nx,1:1:nt,D1,'VA');
  title('Data')
  subplot(3,1,2);
  wiggle(1:1:nx,1:1:nt,D2,'VA');
  title('Observed')
  
  subplot(3,1,3);
  wiggle(1:1:nx,1:1:nt,D2);
  hold on
end
wiggle(1:1:nx,[1:1:nt],D1);
hold off
    