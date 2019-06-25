currentdir = pwd;
cd pvt/;
%CCFLAGS="-march=core2 -O2 -msse -msse2"
mex nuft_forw.cpp
mex nuft_back.cpp
mex nuft_kern.cpp
mex nuft_forw_gridding.cpp
mex nuft_back_gridding.cpp
mex myerfzparts_c.cpp
mex insidepoly_dblengine.c
cd(currentdir);
addpath('pvt/');
