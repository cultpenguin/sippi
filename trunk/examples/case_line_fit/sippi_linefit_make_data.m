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
rand('seed',1);randn('seed',1);

%% Select reference model
m_ref{1}=-30;
m_ref{2}=2;
m_ref{3}=0; 

%% Setup the forward model in the 'forward' structure
nd=40;
forward.x=linspace(1,20,nd);
forward.forward_function='sippi_forward_linefit';

%% Compute a reference set of observed data

x=linspace(1,20,nd)';
d=sippi_forward(m_ref,forward);
d_obs=d{1};
d_std=10;
d_obs=d_obs+randn(size(d_obs)).*d_std;

data{1}.d_obs=d_obs;
data{1}.d_std=d_std;

save sippi_linefit_data forward data m_ref;

%%
plot(forward.x,data{1}.d_obs,'k*')
hold on
e=errorbar(forward.x,data{1}.d_obs,data{1}.d_std.*ones(size(data{1}.d_obs)),'k')
set(e,'LineStyle','none')
hold off
xlabel('x')
ylabel('y')
box on
grid on
axis([0 21 -70 40])
ppp(8,8,12,2,2)
print_mul(sprintf('sippi_linefit_data_%d',nd));
