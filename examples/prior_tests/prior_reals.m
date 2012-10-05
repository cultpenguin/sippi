% Run all the test for generating a priori realizations

cmap=hot;
cmap=(cmap_linear([1 1 0; 0 1 1; 0 0 0]));
cmap=(cmap_linear([1 1 0; 1 0 0; 0 1 0 ;0 0 0]));
cmap=(cmap_linear([1 0 0; 0 1 0 ; 0 0 1 ;0 0 0]));
%cmap=(cmap_linear([1 1 1;1 0 0; 0 1 0 ;0 0 0]));
%cmap=sippi_colormap;
cmap=(cmap_linear([1 0 0; 0 1 0 ;0 0 0]));

save cmap cmap

clear all,close all
prior_reals_fftma
%%
clear all,close all
prior_reals_fftma_cov

clear all, close all
prior_reals_snesim

clear all, close all
prior_reals_visim


sgems_clean;
visim_clean;