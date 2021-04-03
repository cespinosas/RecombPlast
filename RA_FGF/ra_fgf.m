%Differential equations that models the antagonistic gradients of Retinoic Acid (RA) and fibroblast growth factor (FGF) in the presomitic mesoderm of vertebrates as described in Goldbeter et al 2007.

function dydt = ra_fgf(t,y, vs1, Mf, V0, Vsc, kd1, kd2, kd3, kd4, kd5, ks2, ks3, Ka, Ki, n, m)


dydt = [vs1 - kd1*y(3)*y(1) - kd5*y(1);
        V0 + Vsc*(y(4)^n)/(Ka^n + y(4)^n) - kd3*y(2);
        ks2*y(2) - kd2*y(3);
        ks3*Mf*(Ki^m)/(Ki^m + y(1)^m) - kd4*y(4)];
