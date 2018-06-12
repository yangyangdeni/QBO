function [ u ] = hl1972( ~ )
% ���� 201805
% Holton & Lindzen 1972��QBOģ��
% ���ں��������ȴ�Ķ�Kelvin����G-R��������
% ��λ�ϣ�ʱ����1day,z�ĵ�λ��km�������ñ�׼��λ
z=(17:0.25:35)';%250m�ֱ���
n1=length(z);
L=7200;
u=zeros(n1,L);
%delta(t)=1day,2160����6��
S=4.75*10^(-5);
H=6;%���
rio=exp(-z/H);
a=zeros(n1,1);
for i=1:n1
    if z(i)>=30
        a(i)=1/7;
    else
        a(i)=1/21+(1/7-1/21)*(z(i)-17)/13;
    end
end
c1=30*3600*24;k1=2*pi/40000000;%Kelvin
c2=-30*3600*24;k2=2*pi/10000000;%Mixed
for time=2:L;
    u(:,time)=x_cal(u(:,time-1),z,S,rio,a,c1,c2,k1,k2,time);
end
pcolor(0:1/360:20-1/360,17:0.25:35,u/3600/24);colormap(jet)
shading interp;colorbar
xlabel('Year','Fontsize',20)
ylabel('Z(km)','Fontsize',20)
title('Holton & Lindzen 1972 QBO Model','Fontsize',30)
end

