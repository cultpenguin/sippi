% sippi_linefit_make_data: Make data for simple linefit problem
% 
% rand('seed',1);randn('seed',1);
% 
% grad=2;
% intercept=-30;
% poly2=0;
% 
% nd=40;
% x=linspace(1,20,nd)';
% d_obs=x.^2.*poly2+x.*grad+intercept;
% 
% d_std=10;
% d_obs=d_obs+randn(size(d_obs)).*d_std;
% 
% 
% save sippi_linefit_data x d_obs d_std poly2 grad intercept

%%
clear all;close all;
%rand('seed',1);randn('seed',1);
rng('default');rng(1);

%% Select reference model
m_ref{1}=-30;
m_ref{2}=4;
m_ref{3}=-.15; 

%% Setup the forward model in the 'forward' structure
nd=11;
x=linspace(0,20,nd)';
forward.x=x;
forward.forward_function='sippi_forward_linefit';

%% Compute a reference set of observed data

d=sippi_forward(m_ref,forward);
d_ref=d{1};
d_std=zeros(size(d_ref))+10;
d_noise=randn(size(d_ref)).*d_std;
d_obs=d_ref+d_noise;

%data{1}.d_obs=d_obs;
%data{1}.d_std=d_std;

%save sippi_linefit_data forward data m_ref;
save sippi_linefit_data x d_obs d_std m_ref d_ref

%%
figure(1);clf;
plot(x,d_obs,'k*')
hold on
e=errorbar(x,d_obs,d_std,'k');
set(e,'LineStyle','none')
hold off
xlabel('x')
ylabel('d')
box on
grid on
axis([-1 21 -60 40])
ppp(10,7,12,2,2)
print_mul(sprintf('sippi_linefit_data_%d',nd));
