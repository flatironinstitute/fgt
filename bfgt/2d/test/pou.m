function fout=pou(x,center,bsize)

r1=bsize/2*0.9;
r0=r1*0.4;

fout=erfc(12*(sqrt(sum((x-center).^2,1))-(r0+r1)/2)/(r1-r0))/2;
  
end

