% plot_crosshole_data Plot data set from Arrenaes
%
% For more information on the data see 
%
% Hansen, T.M., Cordua, K. S., Looms. M.C., and Mosegaard, K., accepted 2012. SIPPI SIPPI : A Matlab toolbox for Sampling the solution to Inverse Problems with complex Prior Information: Part 2 - Application to cross hole GPR tomography. Computers & Geosciences. doi:10.1016/j.cageo.2012.10.001.
% and 
% Looms, M. C., Hansen, T. M., Cordua, K. S., Nielsen, L., Jensen, K. H., Binley, A., 2010. Geostatistical inference using crosshole ground-penetrating radar : Geostatistical inference using GPR, Geophysics, 75(6), pp J29-J41. doi:10.1190/1.3496001
%
clear all;close all
D{1}='AM13_data';
D{2}='AM24_data';
D{3}='AM1234_data';

for id=1:length(D);

load(D{id})
nd=length(d_obs);
figure(id*10+1);
h=errorbar(1:1:nd,d_obs,d_std,'k*');
nd=length(d_obs);
hold on
plot(1:1:nd,d_obs,'r*')
hold off
errorbarT(h,1,.1)
xlabel('Data #')
ylabel('d_{obs} (ms)')
title(sprintf('%s',D{id}),'interp','none')
print_mul(sprintf('%s',D{id}))

%%
figure(id*10+2);clf;set_paper('portrait');
v_min=0.12;
v_max=0.18;
for i=1:nd
    
    dis=sqrt(sum((R(i,:)-S(i,:)).^2));
    v_app(i)=dis/d_obs(i);
    
    if (size(S,2)==2)
        p=plot([S(i,1),R(i,1)],[S(i,2),R(i,2)],'k-') ;
    else
        p=plot3([S(i,1),R(i,1)],[S(i,2),R(i,2)],[S(i,3),R(i,3)],'k-') ;
    end
        cmap=colormap;
    
    icol=ceil((v_app(i)-v_min)/(v_max-v_min)*size(cmap,1));
    icol=max([icol 1]);
    set(p,'Color',cmap(icol,:));
    if i==1;hold on;end
      
end
axis image
if (size(S,2)==2)
    % 2D
    axis([-.5 5.5 0 12.5]);
    set(gca,'ydir','revers')
    xlabel('X (m)')
    ylabel('Z (m)')
else
    % 3D
    %axis([-.5 3.9 -.5 3.9 0 12.5]);
    set(gca,'zdir','revers')
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Z (m)')
end
hold off
caxis([v_min v_max]);colorbar
title('Apparent velocity (m/ns)')
print_mul(sprintf('%s_raycoverage',D{id}))




end


