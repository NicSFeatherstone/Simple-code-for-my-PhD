//Code for partial pressure of water in reactor after injection of water
//Nicholas Featherstone - 20/10/2023
clear, close, clc

//=======================================================
//Constants
//=======================================================

//Universal Constants
R = 0.0831447                   //bar.mL/mmol.K
Vmol = 22.4                     //mL/mmol
MH2O = 18                       //mol/g

//Reactor Constants
Conversion = 0                  //Conversion of CO2
V = 500                         //mL
T = 280+273.15                  //K
Ptot = 20                       //bar
PCO2In = 1/4.5*Ptot             //Partial Pressure of CO2
ntot = (Ptot*V)/(R*T)           //mmol
PH2O0 = PCO2In*2*Conversion     //bar, partial pressure of water before addition of liquid water, based on conversion of CO2 assuming formation of hydrocarbons (no CO)

//Flow Constants
RestInml = 112.5                //ml/min
RestIn = RestInml/Vmol          //mmol/min

H2OIng = 1                      //g/min = ml/mol
H2OIn = H2OIng/MH2O*1000        //mmol/min
H2OFlowtime = 1                 //min
TotalH2O = H2OIng*H2OFlowtime   //g

//Pressure.Time at pressure greater than...
Pgreater = 1                    //bar

//=======================================================
//Functions
//=======================================================

function f = dpH2Odtup(t,pH2O)    //DE
    f = (Ptot/ntot)*(H2OIn-(RestIn+H2OIn)*(pH2O/Ptot))
endfunction

function f = dpH2Odtdown(t,pH2O)    //DE for changed flow
    f = (Ptot/ntot)*(0-(RestIn+0)*(pH2O/Ptot))
endfunction

//=======================================================
//Solving
//=======================================================

//-------------------------------------------------------
//Solving for no more flow of water
//-------------------------------------------------------

to = 0                          //time water added: 0 min
tup = linspace(to,H2OFlowtime,1000)

Pupvec = ode(PH2O0, to, tup, dpH2Odtup)

PGreaterUp = []
TimeStart = 0
for i = 1:length(tup)
    if Pupvec(i)< Pgreater
        PGreaterUp(i)=0
    else
        PGreaterUp(i) = Pupvec(i)-Pgreater
    end
    
    if PGreaterUp(i)==0
        TimeStart = tup(i+1)
    end
end

//-------------------------------------------------------
//Solving for no more flow of water
//-------------------------------------------------------

tdown = linspace(tup($),1000,1000)
PH2O1 = Pupvec($)

Pdownvec = ode(PH2O1, tup($), tdown, dpH2Odtdown)

PGreaterDown = []
TimeEnd = 0
for i = 1:length(tdown)
    if Pdownvec(i)< Pgreater
        PGreaterDown(i)=0
    else
        PGreaterDown(i) = Pdownvec(i)-Pgreater
    end

    if PGreaterDown(i)>0
        TimeEnd = tdown(i)
    end
end

//=======================================================
//Output
//=======================================================

scf(1),clf
subplot(1,2,1)
plot(logflag="ln",tup,Pupvec)
plot(logflag="ln",tdown,Pdownvec)
//linex = [tup($)tup($)]
//liney = [0,Pupvec($)+0]
//plot(linex,liney,'color','darkgreen')
xlabel('time (min)')
ylabel('PH2O (bar)')
title('log-ln')

subplot(1,2,2)
plot(tup,Pupvec)
plot(tdown,Pdownvec)
//linex = [tup($)tup($)]
//liney = [0,Pupvec($)+0]
//plot(linex,liney,'color','darkgreen')
xlabel('time (min)')
ylabel('PH2O (bar)')
title('ln-ln')

Pxtup = inttrap(tup,Pupvec)      //trapezoidal integration
Pxtdown = inttrap(tdown,Pdownvec)//of pressure * time
Pxttotal = Pxtup+Pxtdown

Pxtupgreater = inttrap(tup,PGreaterUp)      //trapezoidal integration
Pxtdowngreater = inttrap(tdown,PGreaterDown)//of pressure * time
Pxttotalgreater = Pxtupgreater+Pxtdowngreater

TimeGreater = TimeEnd + TimeStart

mprintf('Rest flow rate %.2f ml/min',RestInml)
mprintf('\nWater flow rate %.2f g/min',H2OIng)
mprintf('\nWater flow time %.2f min',H2OFlowtime)
mprintf('\nTotal water in %.2f min',TotalH2O)
mprintf('\nMax partial pressure of water is %.2f bar',Pupvec($))
//mprintf('\nUp time*pp of water is %.2f bar.min',Pxtup)
//mprintf('\nDown time*pp of water is %.2f bar.min',Pxtdown)
mprintf('\nTotal time*pp of water is %.2f bar.min',Pxttotal)
//mprintf('\nUp time*pp >%.2f bar is %.2f bar.min',Pgreater,Pxtupgreater)
//mprintf('\nDown time*pp >%.2f bar is %.2f bar.min',Pgreater,Pxtdowngreater)
mprintf('\nTotal time*pp >%.2f bar is %.2f bar.min',Pgreater,Pxttotalgreater)
mprintf('\nTotal time >%.2f bar is %.2f min',Pgreater,TimeGreater)
