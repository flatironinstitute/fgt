function plotsol()
a=load('fort.34');
x=a(:,1);f=a(:,2);u=a(:,3);rerr=a(:,4);
% f=abs(f);u=abs(u);
% f(f<eps)=eps;
% u(u<eps)=eps;
% semilogy(x,abs(f),'r.')
% hold on;
% semilogy(x,abs(u),'b.')
% semilogy(x,rerr,'c.')
% legend('abs(f)','abs(uex)','abs(uex-ucomp)')
size(f)
[x,I]=sort(x);
f=f(I);
plot(x,f,'r')