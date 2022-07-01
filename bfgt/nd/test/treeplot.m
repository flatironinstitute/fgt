xx = load("tree.data");

clf
figure(1)
hold on;

nlfbox = size(xx,1);
for i = 1:nlfbox
    plot(xx(i,1:5),xx(i,6:10),'k-','LineWidth',0.3)
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

print -depsc bfgttreeout.eps
