%This program takes two random networks that produce two GAPs and generates their offspring by recombination. Then the program tests their offspring accesibility to the atractor B by mutation and plasticity, changing one of the parameters and the initial conditions respectively.

seed = str2num(argv(){1});
nugefo = 500;
nuper = 250;
k2 = 0.25;
k1 = 5*k2;

original_ics = [50, 50];
par_art = [50, 50];
no_form = 2;

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

file_par = fopen(strcat("S", num2str(seed),"/access_ci_mu.txt"), "r");
mut_par = textscan(file_par, "%f %f");
fclose(file_par);

file_out = fopen(strcat("S", num2str(seed),"/pci_pmu_hci_hmu_has.txt"), "w");

"empieza"
for i = 1:nugefo
    da1 = gen_foc{1}(2*i - 1);
    da2 = gen_foc{2}(2*i - 1);
    
    ma1 = gen_foc{1}(2*i);
    ma2 = gen_foc{2}(2*i);
    
    acidad = mut_par{1}(2*i - 1);
    amudad = mut_par{2}(2*i - 1);
    acimom = mut_par{1}(2*i);
    amumom = mut_par{2}(2*i);

    midpci = (acidad + acimom)/2;
    midpmu = (amudad + amumom)/2;

    contBci = 0;
    contBmu = 0;
    
    contBas = 0;
    
    for h = 1:(2^no_form)
    
        if(mod((h-1),2))
            alfa1 = da1;
        else
            alfa1 = ma1;
        end
        
        if(mod(floor((h-1)/2),2))
            alfa2 = da2;
        else
            alfa2 = ma2;
        end
        
        [t,y] = ode45(@(t,y) toggle(t, y, alfa1, alfa2, 3, 3), [0 20],[y1_0f; y2_0f]);
        [m, n] = size(y);
        disA = sqrt(((y(m,1)-alfa1)*(y(m,1)-alfa1))+(y(m,2)*y(m,2)));
        disB = sqrt(((y(m,2)-alfa2)*(y(m,2)-alfa2))+(y(m,1)*y(m,1)));
        if (disA > disB)
            contBas = contBas + 1;
        end
        
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
    end
    
    contBci = contBci/(2^no_form);
    contBmu = contBmu/(2^no_form);
    fprintf(file_out, '%f %f %f %f %f \n', midpci, midpmu, contBci, contBmu, contBas);
    
end
fclose(file_out);

