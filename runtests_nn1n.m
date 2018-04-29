%clear;
%rng('shuffle')
N = 1e6;
gap = 10; skip = 100; 
nMC = gap*N;
beta = 2.5; gamma = 0.8;
dt = 1/(100*beta);
d = 4;
filename ='temp';

x0 = randn(d,1);
p0 = randn(d,1);

tic;

xs = BDsampler_nn1n(datax, datay, ...
        x0, dt, nMC, filename, gap, skip, beta, gamma);

toc;
     
save 'xs_nn1n.mat' xs
