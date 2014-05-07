% sippi_linefit_make_data: Make data for simple linefit problem

rand('seed',1);randn('seed',1);

grad=2;
intercept=-30;
poly2=0;

nd=40;
x=linspace(1,20,nd)';
d_obs=x.^2.*poly2+x.*grad+intercept;

d_std=10;
d_obs=d_obs+randn(size(d_obs)).*d_std;


save sippi_linefit_data x d_obs d_std poly2 grad intercept