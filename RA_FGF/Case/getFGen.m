%This program obtains the random networks that produce two GAPs from the initial conditions of focal_ic.txt by changing the model equation parameters.

seed = str2num(argv(){1});
nugefo = 1000;

rand("seed", seed);
par_art = [4.97, 1.5, 0.365, 7.1, 1.0, 0.28, 1.0, 1.0, 0.0, 1.0, 1.0, 0.2, 0.25];

vs1 = zeros(nugefo, 1);
Mf = zeros(nugefo, 1);
V0 = zeros(nugefo, 1);
Vsc = zeros(nugefo, 1);
kd1 = zeros(nugefo, 1);
kd2 = zeros(nugefo, 1);
kd3 = zeros(nugefo, 1);
kd4 = zeros(nugefo, 1);
kd5 = zeros(nugefo, 1);
ks2 = zeros(nugefo, 1);
ks3 = zeros(nugefo, 1);
Ka = zeros(nugefo, 1);
Ki = zeros(nugefo, 1);
n = m = 2;

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

file_ic = fopen(strcat("S", num2str(seed),"/focal_ic.txt"), "r");
x = textscan(file_ic, "%f");
fclose(file_ic);
y1_0f = x{1}(1,1);
y2_0f = x{1}(2,1);
y3_0f = x{1}(3,1);
y4_0f = x{1}(4,1);

addpath(genpath("../"));

fileout = fopen(strcat("S", num2str(seed),"/genotiposfocales.txt"), "w");

for i = 1:nugefo
    ya = false;
    while(!ya)    
        
        vs1 = par_art(1)*(rand() + 0.5);
        Mf = par_art(2)*(rand() + 0.5);
        V0 = par_art(3)*(rand() + 0.5);
        Vsc = par_art(4)*(rand() + 0.5);
        kd1 = par_art(5)*(rand() + 0.5);
        kd2 = par_art(6)*(rand() + 0.5);
        kd3 = par_art(7)*(rand() + 0.5);
        kd4 = par_art(8)*(rand() + 0.5);
        kd5 = par_art(9)*(rand() + 0.5);
        ks2 = par_art(10)*(rand() + 0.5);
        ks3 = par_art(11)*(rand() + 0.5);
        Ka = par_art(12)*(rand() + 0.5);
        Ki = par_art(13)*(rand() + 0.5);

        [t,y] = ode45(@(t,y) ra_fgf(t,y, vs1, Mf, V0, Vsc, kd1, kd2, kd3, kd4, kd5, ks2, ks3, Ka, Ki, n, m), [0 20],[y1_0f; y2_0f; y3_0f; y4_0f]);
        [msize, nsize] = size(y);
        
        rho = (y(msize,4)/(y(msize,4)+1)) / (y(msize,1)/(y(msize,1)+1));
    
        if (rho < 1/2)
            for c = 1:20
                y1eps = y4_0f + abs(randgauss()*0.15);
                y4eps = y1_0f + abs(randgauss()*0.5);
                [t,y] = ode45(@(t,y) ra_fgf(t,y, vs1, Mf, V0, Vsc, kd1, kd2, kd3, kd4, kd5, ks2, ks3, Ka, Ki, n, m), [0 20],[y1eps; y2_0f; y3_0f; y4eps]);
                [msize, nsize] = size(y);
        
                rho = (y(msize,4)/(y(msize,4)+1)) / (y(msize,1)/(y(msize,1)+1));
                if (rho > 2)
                    ya = true;
                    break;
                end
            end
        end
    endwhile
    fprintf(fileout,'%f %f %f %f %f %f %f %f %f %f %f %f %f \n', vs1, Mf, V0, Vsc, kd1, kd2, kd3, kd4, kd5, ks2, ks3, Ka, Ki);
    clear vs1 Mf V0 Vsc kd1 kd2 kd3 kd4 kd5 ks2 ks3 Ka Ki t y disA disB ya;
end

fclose(fileout);
"listo_get"
