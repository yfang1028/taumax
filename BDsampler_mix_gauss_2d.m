function xs = BDsampler_mix_gauss_2d(x0, dt, nMC, filename, gap, skip)

x = x0;
d = length(x0);

F = force_mix_gauss_2d(x);

fp = fopen(filename,'w');
fclose(fp);

xs = zeros(d,nMC/skip);
for i = 2:nMC
    x = x + dt*F + randn(d,1)*(2*dt)^0.5;
    F = force_mix_gauss_2d(x);
    %{
    if mod(i,gap) == 0
        iter = i/gap;     
        fp = fopen(filename,'a');    
        fprintf(fp,'%d\t%.4f\t%.4f\n',iter,x(1),x(2));
        fclose(fp);    
    end
    %}
    if mod(i,skip) == 0
        xs(:,i/skip) = x;
    end
end


end

