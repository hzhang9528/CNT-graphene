function y=fun2(S)
X=S(1);Y=S(2);h=S(3);alpha1=S(4);beta1=S(4);gamma1=S(4);

fiding=fopen('name.txt','r');
name=fgetl(fiding);
fclose(fiding);


fidin=fopen(name,'r');
a1=fgetl(fidin);
a2=fgetl(fidin);
a3=fgetl(fidin);
a4=fgetl(fidin);
a5=fgetl(fidin);
a6=fgetl(fidin);
fclose(fidin);

aa1 = eval(a1);
aa2 = eval(a2);
aa3 = eval(a3);
aa4 = eval(a4);
aa5 = eval(a5);
aa6 = eval(a6);

y=[aa1;aa2;aa3;aa4;aa5;aa6];

end
