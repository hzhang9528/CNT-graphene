function [subf,X,Y,alpha1,beta1,gamma1,h]=JoiningCNT2(ssgraph,cc,a,n,m,z)

clear subf
syms f X Y alpha1 beta1 gamma1 h

ax=ssgraph(1,:);
ay=ssgraph(2,:);
bx=cc(1,:);
by=cc(2,:);
bz=cc(3,:);

f=0;
for i=1:(m+n)
    f=f+(sqrt((ax(i)+X-(bx(i)*(cos(alpha1)*cos(gamma1)-cos(beta1)*sin(alpha1)*sin(gamma1))+by(i)*(-cos(beta1)*cos(gamma1)*sin(alpha1)-cos(alpha1)*sin(gamma1))+bz(i)*(sin(alpha1)*sin(beta1))))^2 ...
        +(ay(i)+Y-(bx(i)*(cos(gamma1)*sin(alpha1)+cos(alpha1)*cos(beta1)*sin(gamma1))+by(i)*(cos(alpha1)*cos(beta1)*cos(gamma1)-sin(alpha1)*sin(gamma1))+bz(i)*(-cos(alpha1)*sin(beta1))))^2+ ...
        (bx(i)*sin(beta1)*sin(gamma1)+by(i)*cos(gamma1)*sin(beta1)+bz(i)*cos(beta1)+h)^2)-a)^2;
end

a1=diff(f,'X');
a2=diff(f,'Y');
a3=diff(f,'h');
a4=diff(f,'alpha1');
a5=diff(f,'beta1');
a6=diff(f,'gamma1');
functxt=['function',num2str(n),num2str(m),num2str(z),'.txt'];
delete(functxt);
filew = fopen(functxt,'w');
fprintf(filew,'%s\n%s\n%s\n%s\n',char(a1),char(a2),char(a3),char(a4),char(a5),char(a6));
fclose(filew);

nametxt=['name.txt'];
delete(nametxt);
filew2=fopen(nametxt,'w');
fprintf(filew2,['function',num2str(n),num2str(m),num2str(z),'.txt']);
fclose(filew2);
% simplify(a1)
S0=[0.001 0.001 1.95 0.1 0.1 0.1];

% options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e3,'MaxIter',1000,'FunctionTolerance',1e-9);
% options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);
options = optimset('Display','iter','MaxFunEvals',1e4,'MaxIter',1000,'TolFun', 2.0000e-05);
[XYhalpha1beta1gamma1,fval,exitflag,output,jacobian] = fsolve(@fun2, S0,options);
X=XYhalpha1beta1gamma1(1);
Y=XYhalpha1beta1gamma1(2);
h=XYhalpha1beta1gamma1(3);
alpha1=XYhalpha1beta1gamma1(4);
beta1=XYhalpha1beta1gamma1(5);
gamma1=XYhalpha1beta1gamma1(6);
subf=subs(f);

end



