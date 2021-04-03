%This program changes one of the parameters and the initial conditions of the random networks that produce two GAPs, to test if they can access the alternative atractor B by mutation and plasticity, respectively.

seed = str2num(argv(){1});
nugefo = 1000;
nuper = 600;
k2 = 0.25;
k1 = 5*k2;

original_ics = [0.1, 0.1, 0.1, 0.2];
par_art = [4.97, 1.5, 0.365, 7.1, 1.0, 0.28, 1.0, 1.0, 0.0, 1.0, 1.0, 0.2, 0.25];



rand("seed", 1)
function retval = randgauss()
    x = rand();
    y = rand();
    if(x == 0)
        x = rand();
    end
    r = sqrt(-2.0*log(x));
    phi = 2.0*3.141592654*y;
    retval = r*cos(phi);
endfunction

function retval = randgauss_msd(mean, stdev)
    res = randgauss();
    res = res*stdev;
    retval = res+mean;
endfunction

file_ic = fopen(strcat("S", num2str(seed),"/focal_ic.txt"), "r");
x = textscan(file_ic, "%f");
fclose(file_ic);

addpath(genpath("../"));
y1_0f = x{1}(1,1);
y2_0f = x{1}(2,1);
y3_0f = x{1}(3,1);
y4_0f = x{1}(4,1);

file_gen =  fopen(strcat("S", num2str(seed),"/genotiposfocales.txt"), "r");
gen_foc = textscan(file_gen, "%f %f %f %f %f %f %f %f %f %f %f %f %f");
fclose(file_gen);

file_out = fopen(strcat("S", num2str(seed),"/access_ci_mu.txt"), "w");

for i = 1:nugefo
    vs1 = gen_foc{1}(i);
    Mf = gen_foc{2}(i);
    V0 = gen_foc{3}(i);
    Vsc = gen_foc{4}(i);
    kd1 = gen_foc{5}(i);
    kd2 = gen_foc{6}(i);
    kd3 = gen_foc{7}(i);
    kd4 = gen_foc{8}(i);
    kd5 = gen_foc{9}(i);
    ks2 = gen_foc{10}(i);
    ks3 = gen_foc{11}(i);
    Ka = gen_foc{12}(i);
    Ki = gen_foc{13}(i);
    n = m = 2;
    
    contBci = 0;
    contBmu = 0;
    for j = 1:nuper
        sor = rand();
        if(sor < 0.25)
            y10p = y1_0f + randgauss_msd(0, original_ics(1)*k1);
            y20p = y2_0f;
            y30p = y3_0f;
            y40p = y4_0f;
            if(y10p < 0)
                y10p = 0;
            end
        elseif(sor < 0.5)
            y10p = y1_0f;
            y20p = y2_0f + randgauss_msd(0, original_ics(2)*k1);
            y30p = y3_0f;
            y40p = y4_0f;
            if(y20p < 0)
                y20p = 0;
            end
        elseif(sor < 0.75)
            y10p = y1_0f;
            y20p = y2_0f;
            y30p = y3_0f + randgauss_msd(0, original_ics(3)*k1);
            y40p = y4_0f;
            if(y30p < 0)
                y30p = 0;
            end
        else
            y10p = y1_0f;
            y20p = y2_0f;
            y30p = y3_0f;
            y40p = y4_0f + randgauss_msd(0, original_ics(4)*k1);
            if(y40p < 0)
                y40p = 0;
            end
        endif 
        
        [t,y] = ode45(@(t,y) ra_fgf(t,y, vs1, Mf, V0, Vsc, kd1, kd2, kd3, kd4, kd5, ks2, ks3, Ka, Ki, n, m), [0 20],[y10p; y20p; y30p; y40p]);
        [msize, nsize] = size(y);
        
        rho = (y(msize,4)/(y(msize,4)+1)) / (y(msize,1)/(y(msize,1)+1));
        
        if (rho > 2)
            contBci = contBci + 1;
        end
    end
    
    for j = 1:nuper
        p1 = vs1;
        p2 = Mf;
        p3 = V0;
        p4 = Vsc;
        p5 = kd1;
        p6 = kd2;
        p7 = kd3;
        p8 = kd4;
        p9 = kd5;
        p10 = ks2;
        p11 = ks3;
        p12 = Ka;
        p13 = Ki;
        
        sor = rand();
        if(sor < 0.25)
            sor2 = rand();
            if(sor2 < 0.33)
                p1 = vs1 + randgauss_msd(0, par_art(1)*k2);
                if(p1 < 0)
                    p1 = 0;
                end
            elseif(sor2 < 0.66)
                p5 = kd1 + randgauss_msd(0, par_art(5)*k2);
                if(p5 < 0)
                    p = 0;
                end
            else
                p9 = kd5;
                if(p9 < 0)
                    p9 = 0;
                end
            endif
            
        elseif(sor < 0.5)
            sor2 = rand();
            if(sor2 < 0.25)
                p3 = V0 + randgauss_msd(0, par_art(3)*k2);
                if(p3 < 0)
                    p3 = 0;
                end
            elseif(sor2 < 0.50)
                p4 = Vsc + randgauss_msd(0, par_art(4)*k2);
                if(p4 < 0)
                    p4 = 0;
                end
            elseif(sor2 < 0.75)
                p7 = kd3 + randgauss_msd(0, par_art(7)*k2);
                if(p7 < 0)
                    p7 = 0;
                end
            else
                p12 = Ka + randgauss_msd(0, par_art(12)*k2);
                if(p12 < 0)
                    p12 = 0.01;
                end
            endif
            
        elseif(sor < 0.75)
            sor2 = rand();
            if(sor2 < 0.5)
                p10 = ks2 + randgauss_msd(0, par_art(10)*k2);
                if(p10 < 0)
                    p10 = 0;
                end
            else
                p6 = kd2 + randgauss_msd(0, par_art(6)*k2);
                if(p6 < 0)
                    p6 = 0;
                end
            endif
        else
            sor2 = rand();
            if(sor2 < 0.25)
                p11 = ks3 + randgauss_msd(0, par_art(11)*k2);
                if(p11 < 0)
                    p11 = 0;
                end
            elseif(sor2 < 0.50)
                p13 = Ki + randgauss_msd(0, par_art(13)*k2);
                if(p13 < 0)
                    p13 = 0.01;
                end
            elseif(sor2 < 0.75)
                p8 = kd4 + randgauss_msd(0, par_art(8)*k2);
                if(p8 < 0)
                    p8 = 0;
                end
            else
                p2 = Mf + randgauss_msd(0, par_art(2)*k2);
                if(p2 < 0)
                    p2 = 0;
                end
            endif
            
        endif
        
        
        [t,y] = ode45(@(t,y) ra_fgf(t,y, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, n, m), [0 20],[y1_0f; y2_0f; y3_0f; y4_0f]);
        [msize, nsize] = size(y);
        
        rho = (y(msize,4)/(y(msize,4)+1)) / (y(msize,1)/(y(msize,1)+1));
        
        if (rho > 2)
            contBmu = contBmu + 1;
        end
    end
    
    fprintf(file_out, '%f %f \n', contBci, contBmu);
end
fclose(file_out);

"listo_access"
