function F = force_nn1n(w, x, y, beta, gamma)

a = w(1);
b = w(2);
c = w(3);
d = w(4);

tanH = tanh(a*x+b);
dtanH = 1. - tanH.^2;
R = y - c*tanH - d;
dLdd = - R;              % R*(-1)
dLdc = dLdd.*tanH;       % R*(-tanH)
dLdb = dLdd.*(c*dtanH);  % R*(-c*dtanH)
dLda = dLdb.*x;          % R*(-c*dtanH*x)

dLda = beta*sum(dLda)+gamma*a;
dLdb = beta*sum(dLdb)+gamma*b;
dLdc = beta*sum(dLdc)+gamma*c;
dLdd = beta*sum(dLdd)+gamma*d;
Uw = [dLda; dLdb; dLdc; dLdd];
F = -Uw;

end