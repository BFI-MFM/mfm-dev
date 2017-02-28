function [pullval] = sde_pull(model_name)
%computes sde pull: mu/sigma -sigma'/2
%input: model name
%output: pull value
LW = 'linewidth'; lw = 1.6;  FS = 'fontsize';
load (model_name);
pullval = mucheb./sigmacheb-diff(sigmacheb)/2.;


figure; plot(pullval,LW, lw);   
xlabel('X', FS,16)
ylabel('pull)',FS,16);