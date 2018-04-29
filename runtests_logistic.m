%clear;
%rng('shuffle')
N = 1e6;
gap = 1; skip = 1;
nMC = gap*N;
dt = 1/20;
d = size(datax,2);
filename = 'temp';

x0 = randn(d,1);
p0 = randn(d,1);

tic;

xs = BDsampler_logistic(datatrain, labeltrain, ...
        x0, dt, nMC, filename, gap, skip, beta, gamma);

toc;
     
save 'xs_logistic.mat' xs
