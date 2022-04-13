function plot1dtree
xsrc=load('sourcetre');
xs=xsrc(:);xs=unique(xs);xs=sort(xs);
%xs=xs(xs>0);
ns=length(xs);
xtrg=load('targettre');
xt=xtrg(:);xt=unique(xt);xt=sort(xt);
%xt=xt(xt>0);
nt=length(xt);
ys=zeros(size(xs));yt=zeros(size(xt));
%ys=xs;yt=xt;
%loglog(xs,ys,'k.');
plot(xs,ys,'k.');
hold on
%loglog(xt,yt,'ro');
plot(xt,yt,'ro');
ylim([-1,1])