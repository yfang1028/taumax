%clear;
%rng('shuffle')
N = 1e6;
gap = 5; skip = 100;
nMC = N*gap; D = 1; dt = 1/50;
d = 2;
filename ='temp';

x0 = [2;0];
p0 = randn(d,1);


tic;

xs = BDsampler_mix_gauss_2d(x0, dt, nMC, filename, gap, skip);

toc;
     
save 'xs_2d.mat' xs
