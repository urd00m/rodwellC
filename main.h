#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

//#define DEBUG 1

// L1490 code just put into a macro
#define L1490 MGO = MGW;  \
                MF = MF1;  \
                MUGA = MUG1;  \
                N = N + 1;  \
                J = J + 1;  \
                JJ = 1; /* year; */  \
                MFA = MF;  \
                TIS = TI;  \
                TP = ((int)TI/1.0)*1.0+TPI;  \
                TZ1 = TP+TZ4;  \
                TZ2 = TZ1+TZ5;  \
                TZ3 = TZ3E;  \
                QBC = QBC1;  \
                MU = MUD;  \
                TAUP = TP+MUGA*.134*RHOW/MUD-TPI;  \
                TPIW = 168.0;  \
                fprintf(fptr,"\t\t\tYEAR\t=\t%d\n",JJ);  \
                fprintf(fptr,"\t\tSTANDBY OR WATER WITHDRAWAL");  \
                fprintf(fptr,"\n");  \
                print_initial_parameter(fptr, MFA, TWB, MUGA, MUD, RHOW, HS, TI); \
                goto L400;

// N2 code put into macro
#define N2 MGO = MGW; \
                MF = MF2; \
                MUGA = MUG2; \
                N = N+1; \
                JJ = 1; \
                MFA = MF; \
                MU = MUD; \
                TIS = TI; \
                TZ1 = TZ2+TZ6; \
                QBC = QBC2; \
                fprintf(fptr,"\t\t\tYEAR\t=\t%d\n",JJ); \
                fprintf(fptr,"\t\tSTANDBY OR WATER WITHDRAWAL"); \
                fprintf(fptr,"\n"); \
                print_initial_parameter(fptr, MFA, TWB, MUGA, MUD, RHOW, HS, TI); \
                goto L400;

// N3 code put into macro
#define N3 MGO = MGW; \
                MUGA = MUGW; \
                MFA = MFS; \
                N = N+1; \
                JJ = 1; \
                MU = MUD; \
                TIS = TI; \
                QBC = QBC3; \
                TZ2 = TZ1+TZS; \
                fprintf(fptr,"\t\t\tYEAR\t=\t%d\n",JJ); \
                fprintf(fptr,"\t\tSTANDBY OR WATER WITHDRAWAL"); \
                fprintf(fptr,"\n"); \
                print_initial_parameter(fptr, MFA, TWB, MUGA, MUD, RHOW, HS, TI); \
                goto L400;

// N4 code put into a macro
#define N4 MGO = MGW; \
                MUGA = MUGS; \
                MFA = MFS; \
                N = N+1; \
                MU = MUD; \
                JJ = JJ+1; \
                TIS = TI; \
                TZ1 = TZ2+TZ6; \
                QBC = QBC4; \
                fprintf(fptr,"\t\t\tYEAR=%d\n",JJ); \
                fprintf(fptr,"S\t\tUMMER WATER WITHDRAWAL"); \
                fprintf(fptr,"\n"); \
                print_initial_parameter(fptr, MFA, TWB, MUGA, MUD, RHOW, HS, TI); \
                goto L400;

// N5 code put into a macro
#define N5 MGO = MGW; \
                MUGA = MUGW; \
                MFA = MFS; \
                N = N+1; \
                MU = MUD; \
                TZ2 = TZ1+TZS; \
                TIS = TI; \
                QBC = QBC5; \
                fprintf(fptr,"\t\t\tYEAR=%d\n",JJ); \
                fprintf(fptr,"W\t\tINTER WATER WITHDRAWAL"); \
                fprintf(fptr,"\n"); \
                print_initial_parameter(fptr, MFA, TWB, MUGA, MUD, RHOW, HS, TI); \
                goto L400;

void print_initial_parameter(FILE *fptr, double MFA, double TWB, double MUGA, double MUD, double RHOW, double HS, double TI);

