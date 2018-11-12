clc;clear;close all;
addpath(genpath('./dde_biftool'));

for i = 1:3
    if i == 1
        bd = [0.01 15];
    else
        bd = [0.01 0.372 1:3:13];
        bd(1) = 0.01;
    end;
    for k = 1:length(bd)-1
        i,bd(k)
        hopf_br = Bif_kappaPbar(i,bd(k),bd(k+1));
        save(['Bif-' num2str(i) '-' num2str(k) '.mat']);
    end;
end;
for i = 2:3
    hopf_br = Bif_tauPM(i);
    save(['Bif-' num2str(k) '-tauPM.mat']);
end;
