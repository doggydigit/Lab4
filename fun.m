N = 10;

d = -2*ones(N,1);
d1= ones(N-1,1);

diag(d,0) + diag(d1,1) + diag(d1,-1)

S = @(phi) sin(phi);