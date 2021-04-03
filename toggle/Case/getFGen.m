%This program obtains the random networks that produce two GAPs from the initial conditions of focal_ic.txt by changing the model equation parameters.

seed = str2num(argv(){1});
nugefo = 1000;

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


rand("seed", 1)
alfa1 = zeros(nugefo, 1);
alfa2 = zeros(nugefo, 1);

file_ic = fopen(strcat("S", num2str(seed),"/focal_ic.txt"), "r");
x = textscan(file_ic, "%f");
fclose(file_ic);

addpath(genpath("../"));
y1_0f = x{1}(1,1);
y2_0f = x{1}(2,1);

fileout = fopen(strcat("S", num2str(seed),"/genotiposfocales.txt"), "w");

for i = 1:nugefo
    ya = false;
    while(!ya)
        alfa1 = 50*rand() + 25;
        alfa2 = 50*rand() + 25;

        [t,y] = ode45(@(t,y) toggle(t, y, alfa1, alfa2, 3, 3), [0 20],[y1_0f; y2_0f]);
        [m, n] = size(y);
        disA = sqrt(((y(m,1)-alfa1)*(y(m,1)-alfa1))+(y(m,2)*y(m,2)));
        disB = sqrt(((y(m,2)-alfa2)*(y(m,2)-alfa2))+(y(m,1)*y(m,1)));
        if (disA < disB)
            for c = 1:20
                y1eps = abs(randgauss()*0.15);
                y2eps = alfa2 + abs(randgauss()*0.5);
                [t,y] = ode45(@(t,y) toggle(t, y, alfa1, alfa2, 3, 3), [0 20],[y1eps; y2eps]);
                [m, n] = size(y);
                disA = sqrt(((y(m,1)-alfa1)*(y(m,1)-alfa1))+(y(m,2)*y(m,2)));
                disB = sqrt(((y(m,2)-alfa2)*(y(m,2)-alfa2))+(y(m,1)*y(m,1)));
                if (disB < disA)
                    ya = true;
                    break;
                end
            end
        end
    endwhile
    fprintf(fileout,'%f %f \n', alfa1, alfa2);
    clear alfa1 alfa2 t y m n disA disB ya;
end

fclose(fileout);
