close all;
clear all;
%addpath(genpath('./mfiles/'));
addpath(genpath('/home/gabriele/Dropbox/Work/Matlab code/Channel Coding Toolbox/'));

if isunix
   rmpath('/home/gabriele/Dropbox/Work/Matlab code/Channel Coding Toolbox/mexWin'); 
elseif ispc 
   rmpath('/home/gabriele/Dropbox/Work/Matlab code/Channel Coding Toolbox/mexLinux'); 
else
    error('No mex files for MAC');   
end