function testpou
  c1=0;
  c2=0.3;
  delta1=0.002;
  delta2=0.2;
  
  fun1 = @(x) exp(-(x-c1).^2/delta1)+exp(-(x-c2).^2/delta2);
  cf = chebfun(fun1,[-1,1],'splitting','on')
  
  x=linspace(-1,1,1000);
  f1=cf(x);
  f2=f1.*pou(x,c1,0.2);
  f3=f1-f2;
  subplot(3,1,1);plot(x,f1)
  subplot(3,1,2);plot(x,f2)
  subplot(3,1,3);plot(x,f3)
end