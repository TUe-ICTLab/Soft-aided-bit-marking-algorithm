close all;
clear all;
addpath(genpath('./'));


if isunix
elseif ispc 
   rmpath('./mexLinux'); 
else
    error('No mex files for MAC');   
end