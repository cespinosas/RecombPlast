%This program takes two random networks that produce two GAPs and generates their offspring by recombination. Then the program tests their offspring accesibility to the atractor B by mutation and plasticity, changing one of the parameters and the initial conditions respectively.

seed = str2num(argv(){1});
nugefo = 500;
nuper = 250; 
k2 = 0.25;
k1 = 5*k2;

original_ics = [0.1, 0.1, 0.1, 0.2];
par_art = [4.97, 1.5, 0.365, 7.1, 1.0, 0.28, 1.0, 1.0, 0.0, 1.0, 1.0, 0.2, 0.25];
no_form = 4;

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

file_par = fopen(strcat("S", num2str(seed),"/access_ci_mu.txt"), "r");
mut_par = textscan(file_par, "%f %f");
fclose(file_par);

file_out = fopen(strcat("S", num2str(seed),"/pci_pmu_hci_hmu_has.txt"), "w");

for i = 1:nugefo
    dvs1 = gen_foc{1}(2*i - 1);
    dMf = gen_foc{2}(2*i - 1);
    dV0 = gen_foc{3}(2*i - 1);
    dVsc = gen_foc{4}(2*i - 1);
    dkd1 = gen_foc{5}(2*i - 1);
    dkd2 = gen_foc{6}(2*i - 1);
    dkd3 = gen_foc{7}(2*i - 1);
    dkd4 = gen_foc{8}(2*i - 1);
    dkd5 = gen_foc{9}(2*i - 1);
    dks2 = gen_foc{10}(2*i - 1);
    dks3 = gen_foc{11}(2*i - 1);
    dKa = gen_foc{12}(2*i - 1);
    dKi = gen_foc{13}(2*i - 1);
    
    mvs1 = gen_foc{1}(2*i);
    mMf = gen_foc{2}(2*i);
    mV0 = gen_foc{3}(2*i);
    mVsc = gen_foc{4}(2*i);
    mkd1 = gen_foc{5}(2*i);
    mkd2 = gen_foc{6}(2*i);
    mkd3 = gen_foc{7}(2*i);
    mkd4 = gen_foc{8}(2*i);
    mkd5 = gen_foc{9}(2*i);
    mks2 = gen_foc{10}(2*i);
    mks3 = gen_foc{11}(2*i);
    mKa = gen_foc{12}(2*i);
    mKi = gen_foc{13}(2*i);
    n = m = 2;
    
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
            p1 = dvs1;
            p5 = dkd1;
            p9 = dkd5;
        else
            p1 = mvs1;
            p5 = mkd1;
            p9 = mkd5;
        end
        
        if(mod(floor((h-1)/2),2))
            p3 = dV0;
            p4 = dVsc;
            p7 = dkd3;
            p12 = dKa;
        else
            p3 = mV0;
            p4 = mVsc;
            p7 = mkd3;
            p12 = mKa;
        end
        
        if(mod(floor((h-1)/4),2))
            p10 = dks2;
            p6 = dkd2;
        else
            p10 = mks2;
            p6 = mkd2;
        end
        
        if(mod(floor((h-1)/8),2))
            p11 = dks3;
            p13 = dKi;
            p8 = dkd4;
            p2 = dMf;
        else
            p11 = mks3;
            p13 = mKi;
            p8 = mkd4;
            p2 = mMf;
        end
        
        
        [t,y] = ode45(@(t,y) ra_fgf(t,y, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, n, m), [0 20],[y1_0f; y2_0f; y3_0f; y4_0f]);
        [msize, nsize] = size(y);
        
        rho = (y(msize,4)/(y(msize,4)+1)) / (y(msize,1)/(y(msize,1)+1));
        if (rho > 2)
            contBas = contBas + 1;
        end
        
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
            end 
            
            [t,y] = ode45(@(t,y) ra_fgf(t,y, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, n, m), [0 20],[y10p; y20p; y30p; y40p]);
            [msize, nsize] = size(y);
        
            rho = (y(msize,4)/(y(msize,4)+1)) / (y(msize,1)/(y(msize,1)+1));
        
            if (rho > 2)
                contBci = contBci + 1;
            end
        end
        
        
        for j = 1:nuper
            pp1 = p1;
            pp2 = p2;
            pp3 = p3;
            pp4 = p4;
            pp5 = p5;
            pp6 = p6;
            pp7 = p7;
            pp8 = p8;
            pp9 = p9;
            pp10 = p10;
            pp11 = p11;
            pp12 = p12;
            pp13 = p13;
            
            sor = rand();
            if(sor < 0.25)
                sor2 = rand();
                if(sor2 < 0.33)
                    pp1 = p1 + randgauss_msd(0, par_art(1)*k2);
                    if(pp1 < 0)
                        pp1 = 0;
                    end
                elseif(sor2 < 0.66)
                    pp5 = p5 + randgauss_msd(0, par_art(5)*k2);
                    if(pp5 < 0)
                        pp5 = 0;
                    end
                else
                    pp9 = p9;
                    if(pp9 < 0)
                        pp9 = 0;
                    end
                endif
                
            elseif(sor < 0.5)
                sor2 = rand();
                if(sor2 < 0.25)
                    pp3 = p3 + randgauss_msd(0, par_art(3)*k2);
                    if(pp3 < 0)
                        pp3 = 0;
                    end
                elseif(sor2 < 0.50)
                    pp4 = p4 + randgauss_msd(0, par_art(4)*k2);
                    if(pp4 < 0)
                        pp4 = 0;
                    end
                elseif(sor2 < 0.75)
                    pp7 = p7 + randgauss_msd(0, par_art(7)*k2);
                    if(pp7 < 0)
                        pp7 = 0;
                    end
                else
                    pp12 = p12 + randgauss_msd(0, par_art(12)*k2);
                    if(pp12 < 0)
                        pp12 = 0.01;
                    end
                endif
                
            elseif(sor < 0.75)
                sor2 = rand();
                if(sor2 < 0.5)
                    pp10 = p10 + randgauss_msd(0, par_art(10)*k2);
                    if(pp10 < 0)
                        pp10 = 0;
                    end
                else
                    pp6 = p6 + randgauss_msd(0, par_art(6)*k2);
                    if(pp6 < 0)
                        pp6 = 0;
                    end
                endif
            else
                sor2 = rand();
                if(sor2 < 0.25)
                    pp11 = p11 + randgauss_msd(0, par_art(11)*k2);
                    if(pp11 < 0)
                        pp11 = 0;
                    end
                elseif(sor2 < 0.50)
                    pp13 = p13 + randgauss_msd(0, par_art(13)*k2);
                    if(pp13 < 0)
                        pp13 = 0.01;
                    end
                elseif(sor2 < 0.75)
                    pp8 = p8 + randgauss_msd(0, par_art(8)*k2);
                    if(pp8 < 0)
                        pp8 = 0;
                    end
                else
                    pp2 = p2 + randgauss_msd(0, par_art(2)*k2);
                    if(pp2 < 0)
                        pp2 = 0;
                    end
                endif
                
            endif
            
            [t,y] = ode45(@(t,y) ra_fgf(t,y, pp1, pp2, pp3, pp4, pp5, pp6, pp7, pp8, pp9, pp10, pp11, pp12, pp13, n, m), [0 20],[y1_0f; y2_0f; y3_0f; y4_0f]);
            [msize, nsize] = size(y);
            
            rho = (y(msize,4)/(y(msize,4)+1)) / (y(msize,1)/(y(msize,1)+1));
        
            if (rho > 2)
                contBmu = contBmu + 1;
            end
        end
    end
    
    contBci = contBci/(2^no_form);
    contBmu = contBmu/(2^no_form);
    fprintf(file_out, '%f %f %f %f %f \n', midpci, midpmu, contBci, contBmu, contBas);
    
end
fclose(file_out);

"listo_fam"
