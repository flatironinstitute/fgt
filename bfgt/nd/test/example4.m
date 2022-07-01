  xg=load('fort.38');
  yg=load('fort.39');
  rhs=load('fort.40');
  uex=load('fort.41');
  ucomp=load('fort.42');
  uerror=load('fort.43');





  nx=2000;

  xg=reshape(xg,nx,nx);
  yg=reshape(yg,nx,nx);
  rhs=reshape(rhs,nx,nx);
  uex=reshape(uex,nx,nx);
  uerror=reshape(uerror,nx,nx);
%
  fig1=figure(1)
  mesh(xg,yg,rhs)
  fmax=max(rhs);
  fmin=min(rhs);
  
%  axis([-0.5 0.5 -0.5 0.5 0 5])
%  caxis([fmin fmax])

  
%   imagesc(xg,yg,fg,[0 fmax]);      
  colormap(jet)


%
%  contourf(xg,yg,fg,'LineStyle','none');
%  plot(x(1),x(2),'r.');
%  pcolor(xg,yg,fg);
%  imagesc(xg,yg,fg);
%  shading flat;
%  colormap(jet(256));
%  colormap(jet);
%  colorbar;
%  caxis([-5e-4,2e-3]);
%  axis equal;
  xlim([-0.5,0.5]);
  ylim([-0.5,0.5]);
  print(fig1,'bfgtrhs','-depsc');

  fig2=figure(2)
  mesh(xg,yg,uex);
  colormap(jet)
  xlim([-0.5,0.5]);
  ylim([-0.5,0.5]);
  print(fig2,'bfgtuex','-depsc');

  fig3=figure(3)
  mesh(xg,yg,uerror);
  colormap(jet)
  xlim([-0.5,0.5]);
  ylim([-0.5,0.5]);
  print(fig3,'bfgtuerror','-depsc');
  
