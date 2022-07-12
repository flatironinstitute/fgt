  xg=load('fort.38');
  yg=load('fort.39');
  uerror=load('fort.44');
  nx=2000;

  xg=reshape(xg,nx,nx);
  yg=reshape(yg,nx,nx);
  uerror=reshape(uerror,nx,nx);
%
  figure
  mesh(xg,yg,uerror);
  colormap(jet)
  xlim([-0.5,0.5]);
  ylim([-0.5,0.5]);
  print -depsc bfgtuerrorb.eps
  
