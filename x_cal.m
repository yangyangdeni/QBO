function [ u2 ] = x_cal( u1,z,S,rio,a,c1,c2,k1,k2,time)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
ut=u1;F1=u1;F2=u1;F=u1;g1=u1;g2=u1;%开好空间，加速
wsa=2*pi/180;%usa在里面再算
b=4*pi/6400000;%beta
g=9.8*(3600*24)^2;
l=length(u1);
errors=172.8;%为避免无穷，不让风速与两个波速相差小于0.002m/s
for i=1:l
   if abs(u1(i)-c1)<errors
       u1(i)=c1+erors*sign(u1(i)-c1);
   elseif abs(u1(i)-c2)<errors
       u1(i)=c2+errors*sign(u1(i)-c2);
   end
   g1(i)=0.5*sqrt(g*S)*a(i)/k1/(c1-u1(i))^2;
   %g2特别容易忘了括号里的1也有单位
   g2(i)=0.5*sqrt(g*S)*a(i)/(c2-u1(i))^3*b/(k2^3)*(1/3600/24-k2^2/b*(u1(i)-c2));
end
F(1)=0;
for i=2:l
    F1(i)=-2*1000*trapz(z(1:i),g1(1:i));
    F2(i)=-2*1000*trapz(z(1:i),g2(1:i));
    F(i)=(4*10^(-3)*rio(1)*exp(F1(i))-4*10^(-3)*rio(1)*exp(F2(i)))*(3600*24)^2;
end
for i=2:l-1
    if z(i)<=28
        G=0;
    else
        %usa=2*(z(i)-28)*sin(wsa*(time-1))*3600*24;
        %G=wsa*usa;
        G=0;
    end
    ut(i)=-1/rio(i)*(F(i+1)-F(i-1))/500+0.3*3600*24*...
        (u1(i+1)-2*u1(i)+u1(i-1))/(250^2)+G;
end
u2=u1+ut;
u2(1)=0;
%u2(end)=2*7*sin(wsa*time)*3600*24;
u2(end)=u2(end-1);
end

