%Differential equations that models the synthetic network developed by Gardner et al 2000.

function dydt = toggle(t,y, alfa1, alfa2, beta, gama)


dydt = [alfa1/(1+y(2)^beta) - y(1); alfa2/(1+y(1)^gama) - y(2)];
