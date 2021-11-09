function A = Riemann(x,y)
for k=1:length(x)-1
    Ad(k) = y(k+1)*(x(k+1)-x(k));
end
A = 2*sum(Ad);
end