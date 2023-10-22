clear, close, clc
//ASF Distrubtion Code
//Nicholas Featherstone 2022

//k = chain length
//a = chain growth probability, note that it is actually 1-a as a is actually the amount of unreacted material to the length k, in the equation.

//generation of wt% of different species at different chain growths

function f=asf(k,a)
    //$w_a(k)=a^2k(1-a)^k^-^1$
    f = (a^2)*k*(1-a)^(k-1)
endfunction

for i = 1:100   //chain growth probability
    a(i)=(i/100)
    for j = 1:1000  //chain length
        k(j)=j
        asfout(j,i) = asf(k(j),a(i))'
    end
end

//matrix then reads
//a-> 0.1, 0.2, ...
//k=1  X    Y
//k=2  Z 
//...

//csvWrite(asfout,'C:\Users\nicho\OneDrive\Desktop\asf all data')

//----------------------------------------------------------
//sorting data by phase at 25 deg C

asfphase(1,:)=asfout(1,:)+asfout(2,:)+asfout(3,:)+asfout(4,:)//gas
asfphase(2,:)=asfout(5,:)//liquid
asfphase(3,:)=asfout(17,:)//wax
for i=1:12
    asfphase(2,:)=asfphase(2,:)+asfout(i+5,:)//liquid
end
for i=1:983
    asfphase(3,:)=asfphase(3,:)+asfout(i+17,:)//wax
end

scf(1),clf()
subplot(1,3,1)
plot(1-a,asfphase')
xlabel('a')
ylabel('wt%')
legend('Gases','Liquids','Wax',[0.1,0.6])
title('phase at 25C 1 bar')

//----------------------------------------------------------
//sorting data by species group

asflength(1,:)=asfout(1,:)
asflength(2,:)=asfout(2,:)+asfout(3,:)+asfout(4,:)
asflength(3,:)=asfout(5,:)
asflength(4,:)=asfout(18,:)
for i=1:12
    asflength(3,:)=asflength(3,:)+asfout(i+5,:)
end
for i=1:982
    asflength(4,:)=asflength(4,:)+asfout(i+17,:)
end

subplot(1,3,2)
plot(1-a,asflength')
xlabel('a')
ylabel('wt%')
legend('C1','C2-C4','C5-C17','C18+',[0.4,0.9])
title('chain length')
csvWrite(a,'C:\Users\nicho\OneDrive\Desktop\a.csv')
csvWrite(asflength','C:\Users\nicho\OneDrive\Desktop\asf by length.csv')

//----------------------------------------------------------
//data for a single chain growth probability

//enter the chain growth you want data for:
chaingrowth = 0.8
//do not enter 1

subplot(1,3,3)
pointtofind=(1-chaingrowth)*100
toplot=log(asfout(:,pointtofind)'./k')
plot(k([1:20],1)',toplot(1,[1:20]),'.')
xlabel('chain length')
ylabel('log(W/N)')
title('species abundance: fixed a = '+string(chaingrowth))

