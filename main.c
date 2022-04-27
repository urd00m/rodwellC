//#define DEBUG 1

int main(void) {
    int i,J,JJ;
    double TZ3,MGO,QBC,MF,TZ4,QBC1,MUG1,MF1,TZ3E;
    double TZ5,MUG2,QBC2,MF2,TZ6,QBC3,QBC4,QBC5, DEPTH;
    double AL,ALPHAI,BO,CPA,CPI,CPW;
    double DT,EIT,E,FI,GAM,H,HA,HB,HI,HS,HBN,HSN,HSO,KI,MU,MUD,MWG;
    double MFS,MFW,MUGS,MUGW;
    double MGW,N,OMEGA,PI,PL,PM,PLT,PMT,PRWT,QS,QT,QTT,QIT,RA,RHOIS,RHOIM;
    double RHOW,RO,TAUP,TI,TIS,TP,TPI,TPIW,TZ1,TZ2,TF,TICE,TWB,TA,TS,TW,ZS,TZS;
    double D,MFA,MW,MWO,HWB,MWGA,LE,AB,HW,AS,VW,AI,VA;
    double MUGA,MUG,ZP,RHOI,ZPP,ZPS,ASP,MUL,DELH,HP,DP,HWBP;
    double TWP,MWP,HF,VWP,DF,EP,PMP,PLP,VAP,AIP;
    double Q,QI,QB,TAU,RHOA,TAP,FB,FBP,B,BP,BZ,TPW;
    double EI,ESR,PRW,EKT,EK,PMG,PLG,EF,EFI,QITI;

    FILE *fptr;
     
    printf("Hello RASC-AL World!\n");
    char path[200];

    getcwd(path, 200);
    printf("Current working directory: %s\n", path);


 /*
    if ((fptr = fopen("input.txt","r")) == NULL){
        printf("Error: opening file for read\n");

        // Program exits if the file pointer returns NULL.
        return(1);
    }

 
//if you want to update this later to rapidly change the input data without recompiling,
//list the input variables here per the example below
 
    fscanf(fptr,"%lf", &TZ3); //hrs


    fclose(fptr);
*/
    
    // ****** PUT YOUR INPUTS HERE ******
    TZ3=216.0;
    MGO=12186.0;
    QBC=68242.84;
    MF=52073.0;
    TZ4=2064.0;
    QBC1=68242.84;
    MUG1=10.0;
    //MUG2=MUG1;
    MF1=52073.0;
    TZ3E=88000.0;
    TZ5=96.0;
    QBC2=0.0;
    MF2=52073.0;
    TZ6=5088.0;
    //MUG2=1000.0;
    QBC3=68242.84;
    QBC4=68242.84;
    QBC5=68242.84;
    DEPTH=160.0;
    MFS=52073.0;
    MFW=52073.0;
    MUGS=10.0;
    MUGW=10.0;
    TICE=-80.0;
    TWB=68.0;
    

    AL = 0.30; // Firn loss parameter 40
    ALPHAI = .0446; // ft2/hr 41
    BO = 1.1; //
    CPA = 0.199; // BTU /lb-F, Cp air 43
    CPI = 0.5; // Cp ice 44
    CPW = 1.0; // Cp water
    
    DT = 8.333001E-03; // hrs (30 secs) 47
    EIT = 0.0;
    E = 0.0;
    FI = 0.90;
    GAM = 1.0;
    H = 10.0;
    HA = 0.725;
    HB = 60.0;
    HI = 0.725;
    HS = 32.5; // BTU/hr-ft2-F 56
    HBN = 24.0;
    HSN = 32.5;
    HSO = 32.5;
    J = 1;
    KI = 1.28; //BTU/hr-ft-F, ice/firn conductivity 61
    MU = 0.0;
    MUD = 7549.5;
    MWG = 0.0; // gallons, bulb water volume in gallons
    MFS=7549.5;
    MFW=7549.5; //5033.0;
    
    MGW = 1106533.0;
    N = 1;
    OMEGA = 5.399;
    PI = 3.141593;
    PL = 0.0;
    PM = 0.0;
    PLT = 0.0;
    PMT = 0.0;
    PRWT = 0.0;
    
    QS = 0.0;
    QT = 0.0;
    QTT = 0.0;
    QIT = 0.0;
    RA = 1.5; //ft, drill radius
    RHOIS = 45.0; //lbm/ft3, start close-off density of firn
    RHOIM = 57.54; // !lbm/ft3, max firn density
    RHOW = 62.6; // ! lbm/ft3, water density
    RO = RA; // ! ft
    
// TIME PARAMETERS
    TAUP = 0.0;
    TI = 0.0;
    TIS = 0.0;
    TP = 8.0;
    TPI = 8.0;
    TPIW = 8.0;
    TZ1 = 8760.0; // ! 8760 days is one year
    TZ2 = 240.0;
    TZS=TZ1-TZ6;
    
    //TZS = TZ1 - TZ6; // ! Summer duration (days)
// TEMPERATURES
    TF = 32.0;
    TA=TICE;
    TS=TICE;
//124.0;
    TW=TWB;
    
    ZS = pow(((RHOIS - 20.18)/2.4996),(1/0.45));

    D = 2.82843*RO;// !ft, diameter of bulb
    MFA = MF ;
    MW = PI * RA * RA * H * RHOW ; //!lbm, water mass
    MWO = MW ;
    HWB = DEPTH + H; //ft, depth to well bottom
    MWGA = MW / (.134 * RHOW); // ! gallons, convert bulb water mass to volume in gallons
    LE = 144.0 + CPI * (TF - TICE) * OMEGA;
    AB = PI * (D * D)/4.0; // ! ft2, air-water interface area
    HW = H; // ! ft, water depth
    AS = 2.0*PI*D*H/3.0; // ! ft2, water-ice contact area
    VW = PI*D*D*H/8.0; // ! ft3, water volume in bulb
    AI = 2.0 * PI * RA * DEPTH; // ! ft2, air-ice contact area
    VA = PI * RA * RA * DEPTH; // ! ft3, air volume

    if ((fptr = fopen("./output.txt","w")) == NULL){
        printf("Error: opening file for write\n");

        // Program exits if the file pointer returns NULL.
        return(1);
    }
L130:
    fprintf(fptr,"%s","WITHDRAWAL RATE = 100 gal/day\n");
L140:
    fprintf(fptr,"BOILER WATER TEMP DEG F = %f\n",TWB);
L150:
    fprintf(fptr,"BOILER WATER FLOW RATE lbm/hr = %f\n",MF);
L160:
    fprintf(fptr,"CONVECTIVE COEFFICIENT BTU/HR-FT2-F = %f\n",HS);
    fprintf(fptr,"INITIAL DRILL RADIUS FT = %f\n",RA);
    fprintf(fptr,"DEPTH TO TOP OF WATER AT START FT = %f\n", DEPTH);
L180:
    fprintf(fptr,"INITIAL PARABOLIC WATER DIAMETER D FT = %f\n", D);
L191:
    fprintf(fptr,"INITIAL PARABOLIC WATER HEIGHT HW FT = %f\n",HW);
L200:
    fprintf(fptr,"INITIAL WATER TEMP TW DEG F = %f\n",TW);
L201:
    fprintf(fptr,"INITIAL AIR TEMP TA DEG F = %f\n",TA);
L202:
    fprintf(fptr,"INITIAL ICE SURFACE TEMP TS DEG F = %f\n",TS);
L210:
    fprintf(fptr,"AMBIENT ICE TEMP DEG F = %f\n",TICE);
L220:
    fprintf(fptr,"EFFECTIVE LATENT HEAT BTU/LB = %f\n",LE);
L222:
    fprintf(fptr,"\n");
L221:
    fprintf(fptr,"TIME IN HRS, WATER VOL MW GALLONS, ICE AREA AI FT2, AIR VOL VA FT3\n");
    fprintf(fptr,"TIME \tTW \tTA \tTS \tMW \tD \tHW \tHWB \tAI \tVA\n");
L253:
    fprintf(fptr,"%4.2f,\t%4.2f,\t%4.2f,\t%4.2f,\t%4.2f,\t%4.2f,\t%4.2f,\t%4.2f,\t%4.2f,\t%4.2f\n",TI, TW, TA, TS, MWGA, D, HW, HWB, AI, VA);

L260:
    //******************TOP OF LOOP***********************
    
    for(i=1;  i<=11250000; i++) {
        printf("TIME=%f\n",TI);
        if (MWG > MGO) goto L1220; //! bulb water volume .gt. initilaize volume
        if (TI > TZ3) goto L1220; // ! time .gt. formation period
        if (J == 1) goto L280; // ! not sure why we branch here, bulb formation?

L400:
        if (TI < TAUP) { //! not sure what taup is
            MF = 0.0;
            MUG = MUGA;
            MU = MUD;
        } else {
            MF = MFA;
            MUG = 0.0;
            MU = 0.0;
        }
L280:
        ZP = HWB-H/2.0; // ! ft, average bulb depth
        RHOI = 20.18 + 2.4996 * pow(ZP,0.45);
        if(ZP > 394.0) RHOI=RHOIM;
    
    //! compute the change in water depth, h (eq. 7)
L291:
        DELH = 16.0*H*(HS*(TW-TF)-QS)*DT/(RHOI*LE*3.0*(2.0*GAM*H+D));
        HP = H+DELH;
        DP = D+GAM*DELH;
        HWBP = HWB+DELH;
//! assumes full shut-off of water leakage into firn at ZS.
        ZPS = HWB-ZS;
        ASP = 2.0*PI*D*H/3.0;
        
        if(ZPS > H) {
            ASP=0.0;
        } else if (HWB > 25.0){
            ZPP = (ZS+HWB-H)/2.0;
            ASP = 2.0*PI*D*H*(1.0-pow((ZPS/H),1.5))/3.0;
            RHOI = 20.18 + 2.4996 * pow(ZPP,0.45);
        }
L283:
        MUL = AL*ASP*(RHOIS - RHOI); // ! water mass lost to rn
        if(MF == 0.0) goto L284;
        TWB = QBC/(CPW*MF) + TW;
L284:
        TWP = TW+(MF*(TWB-TW)-HS*AS*(TW-TF)*(1.0/CPW+(TW-TF)/LE-QS/(LE*HS))-HA*AB*(TW-TA)/CPW)*DT/MW;

        MWP = MW+(((TW-TF)*HS-QS)*AS/LE-MU-MUL)*DT;
        //printf("MWP=%f,AS=%f,%f,%f,%f\n",MWP,AS,LE,MU,MUL);
        MWG = MWP / (.134 * RHOW);
        VWP = MWP / RHOW;
        HF = sqrt(8.0*VWP*HP/PI)/DP;
        DF = DP*sqrt(HF/HP);
        HW = HF;
        EP = CPW * (TWB - TWP) * MF * DT;
        E = E + EP;
        PMP = MU*DT;
        PM = PM + PMP;
        PLP = MUL*DT;
        PL = PL + PLP;
        AIP = AI+PI*((DP*DP)-(D*D))/4.0 + PI*DP*(HP-HF);
        VAP = VA + PI*((DP*DP)*HP-(DF*DF)*HF)/8.0;
        H = HF;
        D = DF;
        TI = DT + TI;
        Q = HI * (TA - TS);
        QI = Q * DT * AI;
        QT = QT + Q * DT;
        QIT = QIT + QI;
        QB = QT / TI;
        TAU = ALPHAI * TI / (RO * RO);
        RHOA = .4758/ (TA + 460.0);
        TAP = TA+(HA*AB*(TW-TA)+HI*AI*(TS-TA))*DT/(RHOA*VA*CPA);
        printf("TAP=%f\n", TAP);
L418:
        FB = (5.0*(BO*BO*BO))/36.0-BO/4.0+1.0/9.0+(1.0/3.0-BO/2.0)*log(BO)-TAU*(BO-1.0+log(BO));
        FBP = (5.0*(BO*BO))/12.0 - .25-log(BO)/2.0+(1.0/3.0-BO/2.0)/BO-TAU*(1.0+1.0/BO);
        BP = BO - FB /FBP;
        BZ = fabs(BP - BO);
        if(BZ < .0001) goto L425;
        BO = BP;
        goto L418;
L425:
        B = BP;
        BO = BP +.1;
        TS = TICE+QB*RO*(B-1.0)*log(B)/(KI*(B-1.0+log(B)));
        
        
        if(J == 1) goto L1031;
        if(TI > TPW) goto L1130;
L1028:
        if(TI > TP) goto L1131;
        goto L560;
L1031:
        if(TI > TP) goto L1128;
L560:
        HWB = HWBP;
        TW = TWP;
        TA = TAP;
        MW = MWP;
        AS = 2.0*PI*D*H/3.0;
        AB = (PI*(D*D))/4.0;
        AI = AIP;
        VA = VAP;
        if (D > 60.0) goto L1010;
        HS = HSO;
        goto L1040;
L1010:
        HS = HSN;
L1040:
        if(TW < 32.0001) goto L1075;
L1041:
        if(TI > TZ2) goto L1220;
        if(TI > TZ1) goto L1220;
    }
    
    
    goto L1760;
L1075:
    TW = 32.0;
    goto L1041;
L1128:
    fprintf(fptr,"%4.1f,\t%4.2f,\t%4.2f,\t%4.2f,\t%4.1f,\t%4.2f,\t%4.2f,\t%4.2f,\t%4.2f,\t%4.2f\n",TI, TWP, TAP, TS, MWG, D, HW, HWBP, AIP, VAP);
    TP = TP + TPI;
    TPW=TP;
    goto L560;
L1130:
    fprintf(fptr,"%4.1f,\t%4.2f,\t%4.2f,\t%4.2f,\t%4.1f,\t%4.2f,\t%4.2f,\t%4.2f,\t%4.2f,\t%4.2f\n",TI, TWP, TAP, TS, MWG, D, HW, HWBP, AIP, VAP);
    TPW = TPW + TPIW;
    goto L1028;
L1131:
    TP = TP + TPI;
    TAUP = TP +MUGA * .134 * RHOW/MUD - TPI ;
    goto L560;
L1220:
    fprintf(fptr,"%4.1f,\t%4.2f,\t%4.2f,\t%4.2f,\t%4.1f,\t%4.2f,\t%4.2f,\t%4.2f,\t%4.2f,\t%4.2f\n",TI, TWP, TAP, TS, MWG, D, HW, HWBP, AIP, VAP);
L1280:
    fprintf(fptr,"\n");
    EI = E - EIT;
    ESR = EI/(TI-TIS);
    EIT = E;
    PRW = MW-MWO + PM;
    PRWT = PRWT+PRW;
    PLT = PLT+PL;
    PMT = PMT+PM;
    EKT = PRWT*19500.0/E;
    EK = PRW * 19500.0 / EI;
    PMG = PM/(.134*RHOW);
    PM = 0.0;
    PLG = PL/(.134*RHOW);
    PL = 0.0;
    MWO = MW;
    EF = E / 140000.0;
    EFI = EI / 140000.0;
    QITI = QIT - QTT;
    QTT = QIT;
L1340:
    fprintf(fptr,"TOTAL ENERGY INPUT BTU \t\t= \t%e\n",E);
    fprintf(fptr,"SEASONAL ENERGY INPUT BTU \t= \t%e\n",EI);
    fprintf(fptr,"SEASONAL ENERGY INPUT GAL FUEL \t= \t%f\n",EFI);
    fprintf(fptr,"SEASONAL ENERGY RATE BTU/HR \t=\t%f\n",ESR);
    fprintf(fptr,"TOTAL ENERGY INPUT GAL FUEL \t=\t%f\n",EF);
    fprintf(fptr,"AVERAGE LB. WATER PER LB. FUEL \t=\t%f\n",EKT);
    fprintf(fptr,"SEASONAL LB. WATER PER LB. FUEL \t=\t%f\n",EK);
    fprintf(fptr,"ENERGY FROM AIR TO ICE BTU \t=\t%e\n",QIT);
    fprintf(fptr,"SEASONAL ENERGY LOSS, AIR/ICE BTU =\t%e\n",QITI);
    fprintf(fptr,"TOTAL WATER WITHDRAWN GAL \t=\t%f\n",PMT/(.134*RHOW));
    fprintf(fptr,"SEASONAL WATER WITHDRAWN GAL \t=\t%f\n",PMG);
    fprintf(fptr,"TOTAL WATER LOSS GAL \t\t=\t%f\n",PLT/(.134*RHOW));
    fprintf(fptr,"SEASONAL WATER LOSS GAL \t\t=\t%f\n",PLG);
L1430:
    fprintf(fptr,"\n");
    if(N == 1) goto L1490;
    if(N == 2) goto L1204;
    if(N == 3) goto L1540;
    //CCC **** END OF YEAR 1 **** 332
    if(N == 4) goto L1520;
    if(N == 5) goto L1500;
    //CCC **** END OF YEAR 2 **** 335
    if(N == 6) goto L1520;
    if(N == 7) goto L1500;
    //CCC **** END OF YEAR 3 **** 338
    if(N == 8) goto L1520;
    if(N == 9) goto L1500;
    //CCC **** END OF YEAR 4 **** 341
    if(N == 10) goto L1520;
    if(N == 11) goto L1500;
    //CCC **** END OF YEAR 5 **** 344
    if(N == 12) goto L1520;
    if(N == 13) goto L1500;
    //CCC **** END OF YEAR 6 **** 347
    if(N == 14) goto L1520;
    if(N == 15) goto L1500;
    //CCC **** END OF YEAR 7 **** 350
    if(N == 16) goto L1520;
    if(N == 17) goto L1500;
    //CCC **** END OF YEAR 8 **** 353
    if(N == 18) goto L1520;
    if(N == 19) goto L1500;
    //CCC **** END OF YEAR 9 **** 356
    if(N == 20) goto L1520;
    if(N == 21) goto L1500;
    //CCC **** END OF YEAR 1O **** 359
    if(N == 22) goto L1760;
L1490:
    MGO = MGW;
    MF = MF1;
    MUGA = MUG1;
    N = N + 1;
    J = J + 1;
    JJ = 1; // ! year;
    MFA = MF;
    TIS = TI;
    int b=TI/8.0;
    TP = (b)*8.0+TPI;
    TZ1 = TP+TZ4;
    TZ2 = TZ1+TZ5;
    TZ3 = TZ3E;
    QBC = QBC1;
    goto L1210;
L1500:
    MGO = MGW;
    MUGA = MUGW;
    MFA = MFS;
    N = N+1;
    MU = MUD;
    TZ2 = TZ1+TZS;
    TIS = TI;
    QBC = QBC5;
    goto L1553;
L1520:
    MGO = MGW;
    MUGA = MUGS;
    MFA = MFS;
    N = N+1;
    MU = MUD;
    JJ = JJ+1;
    TIS = TI;
    TZ1 = TZ2+TZ6;
    QBC = QBC4;
    goto L1551;
L1540:
    MGO = MGW;
    MUGA = MUGW;
    MFA = MFS;
    N = N+1;
    JJ = 1;
    MU = MUD;
    TIS = TI;
    QBC = QBC3;
    TZ2 = TZ1+TZS;
    goto L1550;
L1204:
    MGO = MGW;
    MF = MF2;
    MUGA = MUG2;
    N = N+1;
    JJ = 1;
    MFA = MF;
    MU = MUD;
    TIS = TI;
    TZ1 = TZ2+TZ6;
    QBC = QBC2;
    goto L1550;
L1210:
    MU = MUD;
    TAUP = TP+MUGA*.134*RHOW/MUD-TPI;
    TPIW = 168.0;
L1550:
    fprintf(fptr,"\t\t\tYEAR\t=\t%d\n",JJ);
    fprintf(fptr,"\t\tSTANDBY OR WATER WITHDRAWAL");
    goto L1555;
L1551:
    fprintf(fptr,"\t\t\tYEAR=%d\n",JJ);
    fprintf(fptr,"S\t\tUMMER WATER WITHDRAWAL");
    goto L1555;
L1553:
    fprintf(fptr,"\t\t\tYEAR=%d\n",JJ);
    fprintf(fptr,"W\t\tINTER WATER WITHDRAWAL");
L1555:
    fprintf(fptr,"\n");
L1580:
    fprintf(fptr,"BOILER WATER FLOW RATE lbm/hr \t=\t%f\n",MFA);
    fprintf(fptr,"BOILER WATER TEMPERATURE DEG F \t=\t%f\n",TWB);
L1610:
    fprintf(fptr,"WATER WITHDRAWAL GAL/DAY \t\t=\t%f\n",MUGA);
    fprintf(fptr,"WITHDRAWAL FLOW RATE GAL/MIN \t=\t%f\n",MUD/(8.04*RHOW));
L1640:
    fprintf(fptr,"CONVECTIVE COEFF AFTER R=30 FT BTU/HR-FT2-F \t=\t%f\n",HS);
L1672:
    fprintf(fptr,"START WITHDRAWAL AT HOUR \t\t=\t%f\n",TI);
    fprintf(fptr,"\n");
    goto L400;
L1760:
    fprintf(fptr,"\n");
L1790:
    fprintf(fptr,"TOTAL ENERGY INPUT BTU = %e\n",E);
L1820:
    fprintf(fptr,"TOTAL ENERGY INPUT GAL FUEL = %f\n",E/140000.);
L1821:
    fprintf(fptr,"TOTAL ENERGY LOSS AIR TO ICE BTU = %e\n",QIT);
L1850:
    fclose(fptr);
    printf("END\n");
    return 0;
}
