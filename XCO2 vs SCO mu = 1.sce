clear,close,clc
//Code for SCO vs XCO2 if mu = 1
//Nicholas Featherstone 2022

Krwgs = 0.02 //this is for 283 C ish
n = 1 //for methane, modify
//n = (0.2*1+0.35*3+0.15*5)/(0.2+0.35+0.15) //=2.85 approx average for iron cat
//n = 1000 //just a high value "like infinity"
F = 3 //Stoich ratio
//F = 1 //RWGS Stoich
//F = 4 //Methane Stoich
mu = 1 //by problem definition

points = 100
Xeq = (1+F-sqrt((1+F)^2-4*F*(1-1/Krwgs)))/(2*(1-1/Krwgs))
printf('Xeq = %.2f',Xeq)
Xeq=linspace(Xeq,Xeq,points-1)'

for i=1:(points-1)
    X(i)=i/points
    //mu(i)=0.1345*100*X(i)-2.3562 //For 250-300 C range of gathered data R^2 = 0.3628
    function f = fmu(S)
        f = (X(i)*S*X(i)*(2-S))/((1-X(i))*(F-X(i)-(2+1/n)*X(i)*(1-S)))*(1/Krwgs)-mu//(i)
        //f = (X(i)*S*X(i)*(2-S))/((1-X(i))*(F-X(i)*(3+1/n)+(2+1/n)*X(i)*S))*(1/Krwgs)-mu
    end
    Sstart = 0 //this is very unstable
    [SCO(i),v(i),info(i)] = fsolve(Sstart,fmu)
end

disp(SCO)

scf(1)
clf()
subplot(1,2,1)
plot(X,SCO,Xeq,SCO)
xlabel('$XCO_2$')
ylabel('$SCO$')
title('XCO2 vs SCO')
subplot(1,2,2)
YCO=X.*SCO
plot(X,YCO,Xeq,YCO)
xlabel('$XCO_2$')
ylabel('$YCO$')
title('XCO2 vs YCO')

//scf(2)
//clf()
//subplot(1,2,1)
//plot(X,v,Xeq,v)
//xlabel('$XCO_2$')
//ylabel('residual')
//title('XCO2 vs residual')
//subplot(1,2,2)
//plot(X,info,Xeq,info)
//xlabel('$XCO_2$')
//ylabel('info')
//title('XCO2 vs info')

//scf(3)
//plot(X,mu)
//xlabel('$XCO_2$')
//ylabel('$\mu$')
//title('XCO2 vs mu')
