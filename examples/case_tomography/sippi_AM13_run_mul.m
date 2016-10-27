% sippi_AM13_run_mul: Invert Arrenaes data using a vairety of prior models
clear all;close all

for i=[1,2,3,4];%1:6;
    for j=[1];
        for k=[0 1];
        for l=[0];
        close all;
        clear prior* options* use* do*;
        use_prior=i;
        use_forward=j;
        use_metropolis=1;
        use_rejection=0;
        n_ite=100000;
        n_reals_out=100;
        doAnneal=k;
        doTempering=l;
        sippi_AM13;
        end
        end
    end
end