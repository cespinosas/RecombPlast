%This program changes one of the parameters and the initial conditions of the random networks that produce two GAPs, to test if they can access the alternative atractor B by mutation and plasticity, respectively.

seed = str2num(argv(){1});
nugefo = 1000;
nuper = 600;
k2 = 0.25;
k1 = 5*k2;

original_ics = [50, 50];
par_art = [50, 50];

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

file_gen =  fopen(strcat("S", num2str(seed),"/genotiposfocales.txt"), "r");
gen_foc = textscan(file_gen, "%f %f");
fclose(file_gen);

file_out = fopen(strcat("S", num2str(seed),"/access_ci_mu.txt"), "w");

for i = 1:nugefo
    alfa1 = gen_foc{1}(i);
    alfa2 = gen_foc{2}(i);
    contBci = 0;
    contBmu = 0;
    for j = 1:nuper
        sor = rand();
        if(sor < 0.50)
            y10p = y1_0f + randgauss_msd(0, original_ics(1)*k1);
            y20p = y2_0f;
        else
            y10p = y1_0f;
            y20p = y2_0f + randgauss_msd(0, original_ics(2)*k1);
        endif  
        
        if(y10p < 0)
            y10p = 0;
        end
        if(y20p < 0)
            y20p = 0;
        end
        [t,y] = ode45(@(t,y) toggle(t, y, alfa1, alfa2, 3, 3), [0 20],[y10p; y20p]);
        [m, n] = size(y);
        disA = sqrt(((y(m,1)-alfa1)*(y(m,1)-alfa1))+(y(m,2)*y(m,2)));
        disB = sqrt(((y(m,2)-alfa2)*(y(m,2)-alfa2))+(y(m,1)*y(m,1)));
        if (disA > disB)
            contBci = contBci + 1;
        end
    end
    
    for j = 1:nuper
        sor = rand();
        if(sor < 0.5)
            a1 = alfa1 + randgauss_msd(0, par_art(1)*k2);
            a2 = alfa2;
        else
            a1 = alfa1;
            a2 = alfa2 + randgauss_msd(0, par_art(2)*k2);
        endif
        
        if(a1 < 0)
            a1 = 0;
        end
        if(a2 < 0)
            a2 = 0;
        end
        [t,y] = ode45(@(t,y) toggle(t, y, a1, a2, 3, 3), [0 20],[y1_0f; y2_0f]);
        [m, n] = size(y);
        disA = sqrt(((y(m,1)-alfa1)*(y(m,1)-alfa1))+(y(m,2)*y(m,2)));
        disB = sqrt(((y(m,2)-alfa2)*(y(m,2)-alfa2))+(y(m,1)*y(m,1)));
        if (disA > disB)
            contBmu = contBmu + 1;
        end
    end
    
    fprintf(file_out, '%f %f \n', contBci, contBmu);
end
fclose(file_out);

