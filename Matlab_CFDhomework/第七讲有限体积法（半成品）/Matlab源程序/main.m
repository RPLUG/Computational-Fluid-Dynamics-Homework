clear;clc;
global gamma;gamma=1.4;dt=0.0001;it=1;
generatemesh;
IC;
  while(it<2000)
      disp(it); 
RK;
 [rho,u,v,p] = ctop(U);
  if(min(min(p))<0||min(min(rho))<0)
      break;
  end
  if(abs(max(max(U(:,:,1)-U0(:,:,1))))<1e-6)
      break;
  end
  it=it+1;
  end
figure (1);
 contour(real(xyCell),imag(xyCell),p,'-m');hold on;
figure (2);
 contour(real(xyCell),imag(xyCell),sqrt(u.^2+v.^2),'-m');hold on;
figure (3);
  plot(xy,'m');hold on;
  plot (xy.','m');hold on;
 plot(xyCell_dummyE,'r');hold on;plot (xyCell_dummyE.','r');hold on;plot(xyCell_dummyW,'r');hold on;plot (xyCell_dummyW.','r');hold on;
 plot(xyCell_dummyN,'r');hold on;plot (xyCell_dummyN.','r');hold on;plot(xyCell_dummyS,'r');hold on;plot (xyCell_dummyS.','r');hold on;   plot(xyCell,'k');hold on;
  plot (xyCell.','k');