function [pullval] = sde_pull(model_name)
%computes sde pull: mu/sigma -sigma'/2
%input: model name
%output: pull value
load (model_name);
LW = 'linewidth'; lw = 1.6;  FS = 'fontsize'; 
gg = linspace(domain_x(1), domain_x(2),200);
pullval = feval(mucheb,gg)./feval(sigmacheb,gg)-feval(diff(sigmacheb)./2,gg);


figure; plot(pullval,LW, lw);   
xlabel('X', FS,16)
ylabel('pull',FS,16);