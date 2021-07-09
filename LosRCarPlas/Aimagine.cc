#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <string>
#include "alea.h"
#include "basics.h"
#include "graphb.h"
#include "graphi.h"
#include "graphc.h"


using std::cout;
using std::cin;
using std::getline;
using std::string;
using std::endl;
using std::ofstream;
using std::ifstream;
using namespace std;

void anova2x2(Basics& bx, int nusim, string factor1, string dir1, string dir2, string idps);
void anova2x2x2(Basics& bx, int nusim);

int main(int argc, char **argv) {
    
    int nusim = 1000;
    ofstream fschi, fsR;
    ifstream fe;
    int seed = 1;
    Alea jacta(seed);
    Basics bas(jacta);
    
    system("mkdir AGraphs");

    anova2x2(bas, nusim, "Recombination", "LPop_wnoise_Asex", "LPop_wnoise_Sex", "");
    anova2x2x2(bas, nusim);
    jacta.close_rng();
    return 0;
}

//compPS_F3X2
//compPS_R2x2
//compPS_R3x2

//abrir ofs y añadir ggplot2, hacer Xgraphs
void anova2x2(Basics& bx, int nusim, string factor1, string dir1, string dir2, string idps) {
    int i,j,k;
    string factor2 = "Founder";
    string dirsal = "AGraphs/"+factor1+"_"+factor2+"_"+idps;
    bx.run_command("mkdir "+dirsal);
    string f1levels[2] = {"without", "with"};
    string f2levels[2] = {"A", "AB"};
    string arch[2][2];
    ofstream fs;
    ifstream fe;
    int *vots;
    vots = new int[nusim];
    
    //para tp
    int ***cubtp;
    bx.create_array(cubtp, 2, 2, nusim);
    arch[0][0] = dir1+"/Graphs/Atp.txt";
    arch[0][1] = dir1+"/Graphs/ABtp.txt";
    arch[1][0] = dir2+"/Graphs/Atp.txt";
    arch[1][1] = dir2+"/Graphs/ABtp.txt";

    for (i=0; i < 2; i++) {
        for (j=0; j < 2; j++) {
            bx.open_ifstream(fe, arch[i][j]);
            for (k=0; k < nusim; k++) {
                fe >> cubtp[i][j][k];
            }
            fe.close();
        }
    }
    bx.open_ofstream(fs, dirsal+"/tp.txt");
    fs << factor1 << "\t" << factor2 << "\tGenerations\n";
    for (i=0; i < 2; i++) {
        for (j=0; j < 2; j++) {
            for (k=0; k < nusim; k++) {
                fs << i << "\t" << j << "\t" << cubtp[i][j][k] << endl;
            }
        }
    }
    fs.close();
    
    
    bx.open_ofstream(fs, dirsal+"/desc_stats_tp.txt");
    fs << "***********************EACH CELL**************************\n";
    for (i=0; i<2; i++) {
        for (j=0; j < 2; j++) {
            fs << "***" << factor1 << "," << f1levels[i] << ":" << factor2 << "," << f2levels[j] << ":\n";
            fs << "Mean: " << bx.get_mean(cubtp[i][j], nusim) << endl;
            fs << "Std Dev: " << bx.get_sample_stddev(cubtp[i][j], nusim) << endl;
            bx.sort(cubtp[i][j], vots, nusim);
            fs << "Median: " << bx.get_midpoint(vots, nusim) << endl;
            fs << "Q1: " << bx.get_q1(vots, nusim) << endl;
            fs << "Q3: " << bx.get_q3(vots, nusim) << endl;
            fs << endl;
        }
    }
    fs.close();
    for (i=0; i<2; i++) {
        for (j=0; j < 2; j++) {
            delete [] cubtp[i][j];
        }
        delete [] cubtp[i];
    }
    delete [] cubtp;
    /// hasta aqui tp
    
    //para t50
    int ***cubt50;
    bx.create_array(cubt50, 2, 2, nusim);
    arch[0][0] = dir1+"/Graphs/At50.txt";
    arch[0][1] = dir1+"/Graphs/ABt50.txt";
    arch[1][0] = dir2+"/Graphs/At50.txt";
    arch[1][1] = dir2+"/Graphs/ABt50.txt";

    for (i=0; i < 2; i++) {
        for (j=0; j < 2; j++) {
            bx.open_ifstream(fe, arch[i][j]);
            for (k=0; k < nusim; k++) {
                fe >> cubt50[i][j][k];
            }
            fe.close();
        }
    }
    bx.open_ofstream(fs, dirsal+"/t50.txt");
    fs << factor1 << "\t" << factor2 << "\tGenerations\n";
    for (i=0; i < 2; i++) {
        for (j=0; j < 2; j++) {
            for (k=0; k < nusim; k++) {
                fs << i << "\t" << j << "\t" << cubt50[i][j][k] << endl;
            }
        }
    }
    fs.close();
    
    
    bx.open_ofstream(fs, dirsal+"/desc_stats_t50.txt");
    fs << "***********************EACH CELL**************************\n";
    for (i=0; i<2; i++) {
        for (j=0; j < 2; j++) {
            fs << "***" << factor1 << "," << f1levels[i] << ":" << factor2 << "," << f2levels[j] << ":\n";
            fs << "Mean: " << bx.get_mean(cubt50[i][j], nusim) << endl;
            fs << "Std Dev: " << bx.get_sample_stddev(cubt50[i][j], nusim) << endl;
            bx.sort(cubt50[i][j], vots, nusim);
            fs << "Median: " << bx.get_midpoint(vots, nusim) << endl;
            fs << "Q1: " << bx.get_q1(vots, nusim) << endl;
            fs << "Q3: " << bx.get_q3(vots, nusim) << endl;
            fs << endl;
        }
    }
    fs.close();
    for (i=0; i<2; i++) {
        for (j=0; j < 2; j++) {
            delete [] cubt50[i][j];
        }
        delete [] cubt50[i];
    }
    delete [] cubt50;
    /// hasta aqui t50
    
    //lo de R
    bx.open_ofstream(fs, dirsal+"/paR.sh");
    fs << "library(ggplot2)\n";
    fs << "library(AICcmodavg)\n";
    
    //para tp
    fs << "a <- read.csv(file =\"" << dirsal << "/tp.txt\", head = TRUE, sep=\"\\t\")\n";
    fs << "sinint <- aov(Generations ~ as.factor(" << factor1 << ") + as.factor(" << factor2 << "), data = a)\n";
    fs << "conint <- aov(Generations ~ as.factor(" << factor1 << ") * as.factor(" << factor2 << "), data = a)\n";
    fs << "model.set <- list(sinint, conint)\n";
    fs << "model.names <- c(\"sinint\", \"conint\")\n";
    fs << "akaike <- aictab(model.set, modnames = model.names)\n";
    fs << "capture.output(akaike, file=\"" << dirsal << "/akaike_tp.txt\", append=TRUE)\n";
    
    fs << "susinint <- summary(sinint)\n";
    fs << "capture.output(susinint, file=\"" << dirsal << "/sinint_tp.txt\", append=TRUE)\n";
    fs << "tusinint <- TukeyHSD(sinint)\n";
    fs << "capture.output(tusinint, file=\"" << dirsal << "/sinint_tp.txt\", append=TRUE)\n";
    
    fs << "suconint <- summary(conint)\n";
    fs << "capture.output(suconint, file=\"" << dirsal << "/conint_tp.txt\", append=TRUE)\n";
    fs << "tuconint <- TukeyHSD(conint)\n";
    fs << "capture.output(tuconint, file=\"" << dirsal << "/conint_tp.txt\", append=TRUE)\n";
    
    
    fs << "a$" << factor1 << " <- factor(a$" << factor1 << ")\n";
    fs << "pdf(file=\"" << dirsal << "/bp_tp.pdf\", width=7.5,height=5.32)\n";
    fs << "par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)\n";
    fs << "ggplot(a, aes(x=" << factor1 << ", y=Generations, fill = as.factor(" << factor2 <<"))) + geom_boxplot(outlier.size=0.5, outlier.alpha=0.1, colour=\"grey20\")+ scale_fill_grey(start=0.8, end=0.4, name=\"" << factor2 << "\", labels = c(expression(italic(\'" << f2levels[0] << "\')),expression(italic(\'" << f2levels[1] << "\')))) + scale_x_discrete(labels=c(\"" << f1levels[0] << "\",\"" << f1levels[1] << "\")) + labs(x=expression(\'" << factor1 //cambiar para popsize
    << "\'), y = expression(\'Generations\'), fill = \"\") + theme_classic() + theme(text = element_text(size=14))\n";
    fs << "dev.off()\n";
    
    //para t50
    fs << "a <- read.csv(file =\"" << dirsal << "/t50.txt\", head = TRUE, sep=\"\\t\")\n";
    fs << "sinint <- aov(Generations ~ as.factor(" << factor1 << ") + as.factor(" << factor2 << "), data = a)\n";
    fs << "conint <- aov(Generations ~ as.factor(" << factor1 << ") * as.factor(" << factor2 << "), data = a)\n";
    fs << "model.set <- list(sinint, conint)\n";
    fs << "model.names <- c(\"sinint\", \"conint\")\n";
    fs << "akaike <- aictab(model.set, modnames = model.names)\n";
    fs << "capture.output(akaike, file=\"" << dirsal << "/akaike_t50.txt\", append=TRUE)\n";
    
    fs << "susinint <- summary(sinint)\n";
    fs << "capture.output(susinint, file=\"" << dirsal << "/sinint_t50.txt\", append=TRUE)\n";
    fs << "tusinint <- TukeyHSD(sinint)\n";
    fs << "capture.output(tusinint, file=\"" << dirsal << "/sinint_t50.txt\", append=TRUE)\n";
    
    fs << "suconint <- summary(conint)\n";
    fs << "capture.output(suconint, file=\"" << dirsal << "/conint_t50.txt\", append=TRUE)\n";
    fs << "tuconint <- TukeyHSD(conint)\n";
    fs << "capture.output(tuconint, file=\"" << dirsal << "/conint_t50.txt\", append=TRUE)\n";
    
    
    fs << "a$" << factor1 << " <- factor(a$" << factor1 << ")\n";
    fs << "pdf(file=\"" << dirsal << "/bp_t50.pdf\", width=7.5,height=5.32)\n";
    fs << "par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)\n";
    fs << "ggplot(a, aes(x=" << factor1 << ", y=Generations, fill = as.factor(" << factor2 <<"))) + geom_boxplot(outlier.size=0.5, outlier.alpha=0.1, colour=\"grey20\")+ scale_fill_grey(start=0.8, end=0.4, name=\"" << factor2 << "\", labels = c(expression(italic(\'" << f2levels[0] << "\')),expression(italic(\'" << f2levels[1] << "\')))) + scale_x_discrete(labels=c(\"" << f1levels[0] << "\",\"" << f1levels[1] << "\")) + labs(x=expression(\'" << factor1 //cambiar para popsize
    << "\'), y = expression(\'Generations\'), fill = \"\") + theme_classic() + theme(text = element_text(size=14))\n";
    fs << "dev.off()\n";
    fs.close();
    bx.run_command("Rscript "+dirsal+"/paR.sh");
}



void anova2x2x2(Basics& bx, int nusim) {
    int i,j,k, l;
    string factor1 = "M";
    string factor2 = "Recombination";
    string factor3 = "Founder";
    string dirsal = "AGraphs/de2"+factor1+"_"+factor2+"_"+factor3;
    bx.run_command("mkdir "+dirsal);
    string f1levels[2] = {"50", "2500"};
    string f2levels[2] = {"without", "with"};
    string f3levels[2] = {"A", "AB"};
    string arch[2][2][2];
    ofstream fs;
    ifstream fe;
    int *vots;
    vots = new int[nusim];
    
    string idf1[2] = {"S", "EL"};
    string idf2[2] = {"Asex", "Sex"};
    string idf3[2] = {"A", "AB"};
    
    
    //para tp
    int ****cubtp;
    cubtp = new int***[2];
    for (i = 0; i < 2; i++) {
        bx.create_array(cubtp[i], 2, 2, nusim);
    }
    for (i=0; i < 2; i++) {
        for (j= 0; j < 2; j++) {
            for (k=0; k < 2; k++) {
                arch[i][j][k] = idf1[i]+"Pop_wnoise_"+idf2[j]+"/Graphs/"+idf3[k]+"tp.txt";
            }
        }
    }
    
    for (i=0; i < 2; i++) {
        for (j= 0; j < 2; j++) {
            for (k=0; k < 2; k++) {
                bx.open_ifstream(fe, arch[i][j][k]);
                for (l=0; l < nusim; l++) {
                    fe >> cubtp[i][j][k][l];
                }
                fe.close();
            }
        }
    }
    bx.open_ofstream(fs, dirsal+"/tp.txt");
    fs << factor1 << "\t" << factor2 << "\t" << factor3 << "\tGenerations\n";
    for (i=0; i < 2; i++) {
        for (j= 0; j < 2; j++) {
            for (k=0; k < 2; k++) {
                for (l=0; l < nusim; l++) {
                    fs << f1levels[i] << "\t" << j << "\t" << k << "\t" << cubtp[i][j][k][l] << endl;
                }
            }
        }
    }
    fs.close();
    
    
    bx.open_ofstream(fs, dirsal+"/desc_stats_tp.txt");
    fs << "***********************EACH CELL**************************\n";
    
    for (i=0; i < 2; i++) {
        for (j= 0; j < 2; j++) {
            for (k=0; k < 2; k++) {
                fs << "***" << factor1 << "," << f1levels[i] << ":" << factor2 << "," << f2levels[j] << ":" << factor3 << "," << f3levels[k] << ":\n";
                fs << "Mean: " << bx.get_mean(cubtp[i][j][k], nusim) << endl;
                fs << "Std Dev: " << bx.get_sample_stddev(cubtp[i][j][k], nusim) << endl;
                bx.sort(cubtp[i][j][k], vots, nusim);
                fs << "Median: " << bx.get_midpoint(vots, nusim) << endl;
                fs << "Q1: " << bx.get_q1(vots, nusim) << endl;
                fs << "Q3: " << bx.get_q3(vots, nusim) << endl;
                fs << endl;
            }
        }
    }
    fs.close();
    for (i=0; i<2; i++) {
        for (j=0; j < 2; j++) {
            for (k=0; k < 2; k++) {
                delete [] cubtp[i][j][k];
            }
            delete [] cubtp[i][j];
        }
        delete [] cubtp[i];
    }
    delete [] cubtp;
    /// hasta aqui tp
    
    //para t50
    
    //para t50
    int ****cubt50;
    cubt50 = new int***[2];
    for (i = 0; i < 2; i++) {
        bx.create_array(cubt50[i], 2, 2, nusim);
    }
    for (i=0; i < 2; i++) {
        for (j= 0; j < 2; j++) {
            for (k=0; k < 2; k++) {
                arch[i][j][k] = idf1[i]+"Pop_wnoise_"+idf2[j]+"/Graphs/"+idf3[k]+"t50.txt";
            }
        }
    }
    
    for (i=0; i < 2; i++) {
        for (j= 0; j < 2; j++) {
            for (k=0; k < 2; k++) {
                bx.open_ifstream(fe, arch[i][j][k]);
                for (l=0; l < nusim; l++) {
                    fe >> cubt50[i][j][k][l];
                }
                fe.close();
            }
        }
    }
    bx.open_ofstream(fs, dirsal+"/t50.txt");
    fs << factor1 << "\t" << factor2 << "\t" << factor3 << "\tGenerations\n";
    for (i=0; i < 2; i++) {
        for (j= 0; j < 2; j++) {
            for (k=0; k < 2; k++) {
                for (l=0; l < nusim; l++) {
                    fs << f1levels[i] << "\t" << j << "\t" << k << "\t" << cubt50[i][j][k][l] << endl;
                }
            }
        }
    }
    fs.close();
    
    
    bx.open_ofstream(fs, dirsal+"/desc_stats_t50.txt");
    fs << "***********************EACH CELL**************************\n";
    
    for (i=0; i < 2; i++) {
        for (j= 0; j < 2; j++) {
            for (k=0; k < 2; k++) {
                fs << "***" << factor1 << "," << f1levels[i] << ":" << factor2 << "," << f2levels[j] << ":" << factor3 << "," << f3levels[k] << ":\n";
                fs << "Mean: " << bx.get_mean(cubt50[i][j][k], nusim) << endl;
                fs << "Std Dev: " << bx.get_sample_stddev(cubt50[i][j][k], nusim) << endl;
                bx.sort(cubt50[i][j][k], vots, nusim);
                fs << "Median: " << bx.get_midpoint(vots, nusim) << endl;
                fs << "Q1: " << bx.get_q1(vots, nusim) << endl;
                fs << "Q3: " << bx.get_q3(vots, nusim) << endl;
                fs << endl;
            }
        }
    }
    fs.close();
    for (i=0; i<2; i++) {
        for (j=0; j < 2; j++) {
            for (k=0; k < 2; k++) {
                delete [] cubt50[i][j][k];
            }
            delete [] cubt50[i][j];
        }
        delete [] cubt50[i];
    }
    delete [] cubt50;
    
    //lo de R
    bx.open_ofstream(fs, dirsal+"/paR.sh");
    fs << "library(ggplot2)\n";
    fs << "library(AICcmodavg)\n";
    
    //para tp
    fs << "a <- read.csv(file =\"" << dirsal << "/tp.txt\", head = TRUE, sep=\"\\t\")\n";
    fs << "sinint <- aov(Generations ~ as.factor(" << factor1 << ") + as.factor(" << factor2 << ") + as.factor(" << factor3 << "), data = a)\n";
    fs << "conint <- aov(Generations ~ as.factor(" << factor1 << ") * as.factor(" << factor2 << ") * as.factor(" << factor3 << "), data = a)\n";
    fs << "model.set <- list(sinint, conint)\n";
    fs << "model.names <- c(\"sinint\", \"conint\")\n";
    fs << "akaike <- aictab(model.set, modnames = model.names)\n";
    fs << "capture.output(akaike, file=\"" << dirsal << "/akaike_tp.txt\", append=TRUE)\n";
    
    fs << "susinint <- summary(sinint)\n";
    fs << "capture.output(susinint, file=\"" << dirsal << "/sinint_tp.txt\", append=TRUE)\n";
    fs << "tusinint <- TukeyHSD(sinint)\n";
    fs << "capture.output(tusinint, file=\"" << dirsal << "/sinint_tp.txt\", append=TRUE)\n";
    
    fs << "suconint <- summary(conint)\n";
    fs << "capture.output(suconint, file=\"" << dirsal << "/conint_tp.txt\", append=TRUE)\n";
    fs << "tuconint <- TukeyHSD(conint)\n";
    fs << "capture.output(tuconint, file=\"" << dirsal << "/conint_tp.txt\", append=TRUE)\n";
    
    
    fs << "a$" << factor2 << " <- factor(a$" << factor2 << ")\n";
    fs << "pdf(file=\"" << dirsal << "/bp_tp.pdf\", width=7.5,height=5.32)\n";
    fs << "par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)\n";
    fs << "ggplot(a, aes(x=" << factor2 << ", y=Generations, fill = as.factor(" << factor3 <<"))) + geom_boxplot(outlier.size=0.5, outlier.alpha=0.1, colour=\"grey20\")+ scale_fill_grey(start=0.8, end=0.4, name=\"" << factor3  << "\", labels = c(expression(italic(\'" << f3levels[0] << "\')),expression(italic(\'" << f3levels[1] << "\')))) + scale_x_discrete(labels=c(\"" << f2levels[0] << "\",\"" << f2levels[1] << "\")) + labs(x=expression(\'" << factor2 //cambiar para M
    << "\'), y = expression(\'Generations\'), fill = \"\") + facet_wrap(~ M, labeller = purrr::partial(label_both, sep = \" = \")) + theme_classic() + theme(text = element_text(size=14))\n";
    fs << "dev.off()\n";
    
    //para t50
    fs << "a <- read.csv(file =\"" << dirsal << "/t50.txt\", head = TRUE, sep=\"\\t\")\n";
    fs << "sinint <- aov(Generations ~ as.factor(" << factor1 << ") + as.factor(" << factor2 << ") + as.factor(" << factor3 << "), data = a)\n";
    fs << "conint <- aov(Generations ~ as.factor(" << factor1 << ") * as.factor(" << factor2 << ") * as.factor(" << factor3 << "), data = a)\n";
    fs << "model.set <- list(sinint, conint)\n";
    fs << "model.names <- c(\"sinint\", \"conint\")\n";
    fs << "akaike <- aictab(model.set, modnames = model.names)\n";
    fs << "capture.output(akaike, file=\"" << dirsal << "/akaike_t50.txt\", append=TRUE)\n";
    
    fs << "susinint <- summary(sinint)\n";
    fs << "capture.output(susinint, file=\"" << dirsal << "/sinint_t50.txt\", append=TRUE)\n";
    fs << "tusinint <- TukeyHSD(sinint)\n";
    fs << "capture.output(tusinint, file=\"" << dirsal << "/sinint_t50.txt\", append=TRUE)\n";
    
    fs << "suconint <- summary(conint)\n";
    fs << "capture.output(suconint, file=\"" << dirsal << "/conint_t50.txt\", append=TRUE)\n";
    fs << "tuconint <- TukeyHSD(conint)\n";
    fs << "capture.output(tuconint, file=\"" << dirsal << "/conint_t50.txt\", append=TRUE)\n";
    
    
    fs << "a$" << factor2 << " <- factor(a$" << factor2 << ")\n";
    fs << "pdf(file=\"" << dirsal << "/bp_t50.pdf\", width=7.5,height=5.32)\n";
    fs << "par(cex.lab=1.5, cex.axis=1.5, las=1, mar=c(5,5,2,2)+0.1)\n";
    fs << "ggplot(a, aes(x=" << factor2 << ", y=Generations, fill = as.factor(" << factor3 <<"))) + geom_boxplot(outlier.size=0.5, outlier.alpha=0.1, colour=\"grey20\")+ scale_fill_grey(start=0.8, end=0.4, name=\"" << factor3 << "\", labels = c(expression(italic(\'" << f3levels[0] << "\')),expression(italic(\'" << f3levels[1] << "\')))) + scale_x_discrete(labels=c(\"" << f2levels[0] << "\",\"" << f2levels[1] << "\")) + labs(x=expression(\'" << factor2  << "\'), y = expression(\'Generations\'), fill = \"\") + facet_wrap(~ M, labeller = purrr::partial(label_both, sep = \" = \")) + theme_classic() + theme(text = element_text(size=14))\n";
    fs << "dev.off()\n";

    //para t50
    fs.close();
    bx.run_command("Rscript "+dirsal+"/paR.sh");
}
