//FTHNIC002 - Nicholas Featherstone
//CHE4045Z
//RWGS Eq Code
close,clear,clc, lines(0,100);

// Rxn --> $CO + H_2O <-> H_2 + CO_2$
Tref = 25+273.15 //K
R=8.314

num = 1
while num<=10

    //====================================================================
    //Feed Information: 
    //    CO    CO2   H2    H2O
    F0 = [0     1     num/2   0]';  //mmol/s
    Ftot = sum(F0);
    y0 = F0/Ftot;       //mol fraction
    theta = y0/y0(2)
    thetaCO = theta(1);
    thetaCO2 = theta(2);
    thetaH2 = theta(3);
    thetaH2O = theta(4);

    //====================================================================
    //Thermodynamic Data
    //        CO        CO2         H2   H2O
    DHf298 = [-110530   -393510     0    -241835]; //Heat of fomation at 298K, J/mol
    DGf298 = [-137160   -394380     0    -228614]; //GibbsEnergy of formation at 298K, J/mol
    DHrxn0 = (DHf298(4) + DHf298(1))-(DHf298(2)+DHf298(3)) // Heat of reaction at 25C, J/mol
    DGrxn0 = (DGf298(4) + DGf298(1))-(DGf298(2)+DGf298(3)) //Delta G of formation of reaction at 25C, J/mol

    //Cp data, J/mol.K
    CpCO  = [30.87 -12.85e-3 2.789e-5  -11.68e-9];
    CpCO2 = [19.80 73.44e-3  -5.602e-5 17.15e-9];
    CpH2  = [27.14 9.274e-3  -1.381e-5 7.645e-9];
    CpH2O = [32.24 1.924e-3  1.055e-5  -3.596e-9];
    CpDrxn = CpCO + CpH2O - CpH2 - CpCO2;

    function CpDT = CpT(T1, T2, Cp)
        Tmultiplier = [(T2-T1) (T2^2-T1^2)/2 (T2^3-T1^3)/3 (T2^4-T1^4)/4]'
        CpDT = Cp*Tmultiplier; 
    endfunction

    function f = Cprxn(T)
        f = CpT(Tref,T,CpDrxn)
    endfunction

    function HrxnT= DHrxnT(T)
        HrxnT = DHrxn0 + Cprxn(T)
    endfunction

    //====================================================================
    //Kinetic Data
    function Gibbs = DrxnG(T)
        Tmultiplier = [(Tref) (Tref^2)/2 (Tref^3)/3 (Tref^4)/4]'
        A = DHrxn0 - CpDrxn*Tmultiplier
        Gibbs = T*((DGrxn0/Tref)+A*(1/T-1/Tref)-CpDrxn(1)*log(T/Tref)-CpDrxn(2)*(T-Tref)/2 - CpDrxn(3)*(T^2-Tref^2)/6 - CpDrxn(4)*(T^3-Tref^3)/12)
    endfunction

    function EquilibriumConstant = Ke(T)
        //$ln(K_e_q^0) = -\Delta G_0/RT$
        EquilibriumConstant = exp(-DrxnG(T)/(R*T))
    endfunction

    //Quadratic Formula
    function f4 = Xe(T)
        a=Ke(T)-1;
        b=-(Ke(T)*(thetaH2+1))
        c=Ke(T)*thetaH2
        f4 = (-b-sqrt(b^2-4*a*c))/(2*a)
    endfunction

    //====================================================================

    Tl = 400; //K
    Th = 1500; //K
    T=linspace(Tl,Th,1000)';
    Tin = 400 + 273.15;

    for i=1:length(T)
        Xeq(num,i) = Xe(T(i));
    end
    num=num+1;
end
disp(Xeq(10,:)')
disp(T-273.15)
scf(1),clf()
plot(T-273.15, Xeq(1,:)',T-273.15, Xeq(2,:)',T-273.15, Xeq(3,:)',T-273.15, Xeq(4,:)',T-273.15, Xeq(5,:)',T-273.15, Xeq(6,:)',T-273.15, Xeq(7,:)',T-273.15, Xeq(8,:)');
legend('H2/CO2=0.5','H2/CO2=1','H2/CO2=1.5','H2/CO2=2','H2/CO2=2.5','H2/CO2=3','H2/CO2=3.5','H2/CO2=4',4)
xlabel('C')
ylabel('X')
