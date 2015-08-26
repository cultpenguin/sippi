% zoepprits_forward: Achenbach/Zoeppritz solution
%
% Call:
%   [rt,angle]=zoepprits_forward(Vp1,Vs1,Rho1,Vp2,Vs2,Rho2,angle);
function [rt,angle]=zoepprits_forward(Vp1,Vs1,Rho1,Vp2,Vs2,Rho2,angle);

angle=pi*angle./180;

x=sin(angle);

a=Rho2/Rho1;
b=Vs1/Vp1;
c=Vp2/Vp1;
d=Vs2/Vp1;

M=zeros(4,4);
n=zeros(4,1);

rt=zeros(length(x),4);

for i=1:length(x)

    M(1,1)=-x(i);
    M(2,1)=sqrt(1-x(i)^2);
    M(3,1)=2*(b^2)*x(i)*sqrt(1-x(i)^2);
    M(4,1)=-(1-2*(b^2)*(x(i)^2));
    M(1,2)=-sqrt(1-(b^2)*(x(i)^2));
    M(2,2)=-b*x(i);
    M(3,2)=b*(1-2*(b^2)*(x(i)^2));
    M(4,2)=2*(b^2)*x(i)*sqrt(1-x(i)^2);
    M(1,3)=c*x(i);
    M(2,3)=sqrt(1-(c^2)*(x(i)^2));
    M(3,3)=2*a*(d^2)*x(i)*sqrt(1-(c^2)*(x(i)^2));
    M(4,3)=a*c*(1-2*(d^2)*(x(i)^2));
    M(1,4)= -sqrt(1-(d^2)*(x(i)^2));
    M(2,4)=d*x(i);
    M(3,4)=-a*d*(1-2*(d^2)*(x(i)^2));
    M(4,4)=2*a*(d^2)*x(i)*sqrt(1-(d^2)*(x(i)^2));


    n(1)=x(i);
    n(2)=sqrt(1-x(i)^2);
    n(3)=2*(b^2)*x(i)*sqrt(1-x(i)^2);
    n(4)=1-2*(b^2)*(x(i)^2);
    rt(i,1:4)=(M\n)';

end
