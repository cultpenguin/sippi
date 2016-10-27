% sippi_AM13_run_mul: Invert Arrenaes data using a vairety of prior models
clear all;close all

for i=[1,2,3,4];%1:6;
    for j=[1];
        for an=[0 1];
        for pt=[0 1];
        close all;
        clear prior* options* use* do*;
        use_prior=i;
        use_forward=j;
        use_metropolis=1;
        use_rejection=0;
        n_ite=200000;
        n_reals_out=100;
        doAnneal=an;
        doTempering=pt;
        sippi_AM13;
        end
        end
    end
end