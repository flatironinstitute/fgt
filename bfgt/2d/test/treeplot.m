xx = load("tree.data");
src = load("src.data");
targ = load("targ.data");

clf
figure(1)
hold on;

nlfbox = size(xx,1);
for i = 1:nlfbox
    plot(xx(i,1:5),xx(i,6:10),'k-','LineWidth',0.3)
end

plot(src(:,1),src(:,2),'r-','LineWidth',5)

if ~isempty(targ)
  plot(targ(:,1),targ(:,2),'b.','MarkerSize',8)
end

boxtype = load('type.data');
for i=1:size(boxtype,1)
  if boxtype(i,3) == -3
    text(boxtype(i,1),boxtype(i,2),int2str(boxtype(i,3)),'FontSize',6, 'Color','r');
  else 
    text(boxtype(i,1),boxtype(i,2),int2str(boxtype(i,3)),'FontSize',6);
  end
end

nver=load('nver.data');
vert=load('vert.data');
% plot out interior polygons
% ncomp=length(nver)
% is=1;
% for icomp=1:ncomp
%     ie=is+nver(icomp)-1;
% %    plot(vert(is:ie,1),vert(is:ie,2),'b-')
%     plot(vert(is:ie,1),vert(is:ie,2),'b-','LineWidth',1.5)
%     is=ie+1;
% end

a=load('edge.data');

% nc=size(a,1)/2;
% v1=a(1:nc,:);
% v2=a(nc+1:end,:);
% 
% for i=1:nc
%     plot([v1(i,1) v2(i,1)], [v1(i,2) v2(i,2)], 'r','LineWidth',1.5)
% end
nchs=load('nch.data');
ncomp=length(nchs);

nc=size(a,1)
a=[a; a(1,:)];
is=1;       
ic=0;
for icomp=1:ncomp
    for i=1:nchs(icomp)
        if i<nchs(icomp)
            in=i+1;
        else
            in=1;
        end
        ie=is+nver(ic+i)-1;
        v=[a(ic+i,:); vert(is:ie,:); a(ic+in,:)];
        plot(v(:,1),v(:,2),'r-','LineWidth',1.5)
        is=ie+1;
        %pause
    end
    ic=ic+nchs(icomp);
end

axis equal
axis off


ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];