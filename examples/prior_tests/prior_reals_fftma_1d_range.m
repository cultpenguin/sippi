% prior_reals_fftma_1d_range
clear all;

% prior_reals_fftma_1d_range
im=1; 
prior{im}.type='FFTMA';
prior{im}.x=[0:.1:20]; % X array 
prior{im}.m0=10;
prior{im}.Va='1 Gau(10)';
prior{im}.fftma_options.constant_C=0;


im=2;
prior{im}.type='gaussian';
prior{im}.name='range_1';
prior{im}.m0=5;
prior{im}.std=5;
prior{im}.norm=80;
prior{im}.prior_master=1; % -- NOTE, set this to the FFT-MA type prior for which this prior type
                          % should desriibe the range
                          
                          
%%
figure(1);set_paper('portrait');
clear L mm
randn('seed',1)
rand('seed',1)
m=sippi_prior(prior);
for i=1:10;
    [m,prior]=sippi_prior(prior);
    mm(:,i)=m{1};
    L{i}=format_variogram(prior{1}.Va,1);
end
plot(prior{1}.x,mm,'-');
legend(L,'Location','EastOutside');
xlabel('x(m)')
ylabel('y');
ppp(8,8,12,2,2)    
print_mul('prior_reals_fftma_1d_range')


%% movie
prior{1}.seq_gibbs.step=.005;
prior{2}.seq_gibbs.step=0.1;
aviobj = avifile('prior_reals_fftma_1d_range.avi')

fig=figure(100);

clear m_real;
for i=1:250;
    [m,prior]=sippi_prior(prior,m);
    m_real(:,i)=m{1};
    sippi_plot_prior(prior,m);
    set(gca,'ylim',[7 13])
    drawnow
        F = getframe(fig);
    aviobj = addframe(aviobj,F);

end
aviobj = close(aviobj);

%%
figure(2);clf;set_paper('portrait');
x=prior{1}.x;
ii=10:10:size(m_real,2);
ii=1:1:150;
waterfall(x,ii,m_real(:,ii)');
ylabel('Iteration number')
xlabel('X')
zlabel('Y')
ppp(12,8,12,2,2)
view([-50 40])
print_mul('prior_reals_fftma_1d_range_waterfall')
%close all;



