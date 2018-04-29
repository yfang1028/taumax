function xs = BDsampler_gauss_1d(x0, dt, nMC, filename, gap, skip)

x = x0;
d = length(x0);

F = -x;

fp = fopen(filename,'w');
fclose(fp);

xs = zeros(d,fix(nMC/skip));
for i = 2:nMC
    x = exp(-dt)*x + randn(d,1)*(1-exp(-2*dt))^0.5;    
    %{
    if mod(i,gap) == 0
        iter = i/gap;
        fp = fopen(filename,'a');   
        fprintf(fp,'%d\t%.4f\n',iter,x);
        fclose(fp);    
    end
    %}
    if mod(i,skip) == 0
        xs(:,fix(i/skip)) = x;
    end
end


end

