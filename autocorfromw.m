function tau = autocorfromw(c, w)

M = length(w)+1;
tau = 1 + 2*w'*c(2:M);