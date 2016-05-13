% plot_traveltime_sr
%
% Call 
%    plot_traveltime(S,R)
%    S: [n,2] : source locattion
%    R: [n,2] : reveiver locattion
%    or (3d)
%    S: [n,3] : source locattion
%    R: [n,3] : reveiver locattion
%
%
% EX:
%  % 2D
%  D=load('AM13_data.mat');
%  plot_traveltime_sr(D.S,D.R);
%  or
%  ant_pos=[D.S,D.R];
%  plot_traveltime_sr(ant_pos);
%
%  % 3D
%  D=load('AM1234_data.mat');
%  plot_traveltime_sr(D.S,D.R);
%
%
function plot_traveltime_sr(S,R);

if nargin==1;
    ant_pos=S;
    S=ant_pos(:,1:2);
    R=ant_pos(:,1:2);
end

if size(S,2)==2;
    plot(S(:,1),S(:,2),'r*');
    hold on
    plot(R(:,1),R(:,2),'go');
    for i=1:size(S,1);
        plot([S(i,1) R(i,1)],[S(i,2) R(i,2)],'k-');
    end
elseif size(S,2)==3;
    for i=1:1:size(S,1);
        plot3([S(i,1) R(i,1)],[S(i,2) R(i,2)],[S(i,3) R(i,3)],'k-','LineWidth',.1);
    end
    hold off
end
drawnow;





