function [subf,X,Y,sita,h]=JoiningCNT(ssgraph,cc,a,n,m,z)

clear subf
syms f X Y sita h

ax=ssgraph(1,:);
ay=ssgraph(2,:);
bx=cc(1,:);
by=cc(2,:);

f=0;
for i=1:(m+n)
    f=f+(sqrt((ax(i)+X-(bx(i)*cos(sita)-by(i)*sin(sita)))^2+(ay(i)+Y-(bx(i)*sin(sita)+by(i)*cos(sita)))^2+h^2)-a)^2;
end

a1=diff(f,'X');
a2=diff(f,'Y');
a3=diff(f,'h');
a4=diff(f,'sita');
functxt=['function',num2str(n),num2str(m),num2str(z),'.txt'];
delete(functxt);
filew = fopen(functxt,'w');
fprintf(filew,'%s\n%s\n%s\n%s\n',char(a1),char(a2),char(a3),char(a4));
fclose(filew);

nametxt=['name.txt'];
delete(nametxt);
filew2=fopen(nametxt,'w');
fprintf(filew2,['function',num2str(n),num2str(m),num2str(z),'.txt']);
fclose(filew2);
% simplify(a1)
S0=[0.001 0.001 1.95 0.1];

% options = optimoptions('fsolve','Display','iter','MaxFunEvals',1e3,'MaxIter',1000,'FunctionTolerance',1e-9);
% options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);
options = optimset('Display','iter','MaxFunEvals',1e4,'MaxIter',1000,'TolFun', 2.0000e-05);
[XYhsita,fval,exitflag,output,jacobian] = fsolve(@fun, S0,options);
X=XYhsita(1);
Y=XYhsita(2);
h=XYhsita(3);
sita=XYhsita(4);
subf=subs(f);
