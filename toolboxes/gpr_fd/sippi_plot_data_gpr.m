% sippi_plot_data_gpr: plot GPR data
%                  overwrites sippi_plot_data plot data in SIPPI
%
% Call.
%    sippi_plot_data_gpr(d,data);
%    sippi_plot_data_gpr(d,data,id_arr);
%
function sippi_plot_data_gpr(d,data,id_arr);

if nargin<3, 
    id_arr=1:length(d);
end

figure_focus(20+1);clf;set_paper('portrait')    
for id=id_arr;
    
  D1(:,id)=d{id};
  
  if nargin>1;
    D2(:,id)=data{id}.d_obs;
  end
 
end
[nt,nx]=size(D1);

showmax=300;
scale=[];
type='image';
plImage=1;

if nargin>1
  subplot(3,1,1);
  title('Data')
  wiggle(1:1:nx,1:1:nt,D1,type,scale,showmax,plImage);
  title('Data (forward response)')
  %xlabel('i_{data}'),ylabel('i_t')

  subplot(3,1,2);
  wiggle(1:1:nx,1:1:nt,D2,type,scale,showmax,plImage);
  %xlabel('i_{data}');
  ylabel('i_t')
  title('Data (observed)')
  
  subplot(3,1,3);
  wiggle(1:1:nx,1:1:nt,D1-D2,type,scale,showmax,plImage);
  %hold on
  %wiggle(1:1:nx,1:1:nt,D1,type,scale,showmax,0);
  %hold off
  xlabel('i_{data}');
  
else
    wiggle(1:1:nx,1:1:nt,D1,type,scale,showmax,plImage);
    title('Data (forward response)')
    xlabel('i_{data}'),ylabel('i_t')

end
hold off
    