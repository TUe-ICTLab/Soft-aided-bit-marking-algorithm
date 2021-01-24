%close all;
clear all;
if isunix
addpath(genpath('/home/gliga/Dropbox/Work/Matlab code/Channel Coding Toolbox'));
%addpath(genpath('/home/gabriele/Dropbox/Work/Matlab code/Channel Coding Toolbox'));

elseif ispc 
   %addpath(genpath('E:\Dropbox\Work\Matlab code\Channel Coding Toolbox'));
   addpath(genpath('E:\Users\Gabriele\Dropbox\Work\Matlab code\Channel Coding Toolbox'));
else 
    error('Unexpected Operative System');
end

%addpath(genpath('./SCC/'));