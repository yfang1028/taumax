function xs = BDsampler_nn1n(datax, datay, ...
    x0, dt, nMC, filename, gap, skip, beta, gamma)

x = x0;
d = length(x0);

F = force_nn1n(x, datax, datay, beta, gamma);

fp = fopen(filename,'w');
fclose(fp);

xs = zeros(d,nMC/skip);
for i = 2:nMC
    x = x + dt*F + randn(d,1)*(2*dt)^0.5;
    F = force_nn1n(x, datax, datay, beta, gamma);
    %{
    if mod(i,gap) == 0
        iter = i/gap;
        ff = x;      
        fp = fopen(filename,'a');    
        fprintf(fp,'%d\t%.4f\t%.4f\t%.4f\t%.4f\n', ...
            iter,ff(1),ff(2),ff(3),ff(4));
        fclose(fp);    
    end
    %}
    if mod(i,skip) == 0
        xs(:,i/skip) = x;
    end
end


end

