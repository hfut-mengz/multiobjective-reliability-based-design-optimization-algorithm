function [fval,gnum,dgnum] = PMA(xval,betat,nr,ns,ng)
[mu,sigma,type]=distribution;
n=length(xval);

if ns>0
nsi=n-ns+1;
 for i=nsi:1:n
     mu(i)=xval(i);
 end
else
end

gnum=0;
dgnum=0;
for i=1:1:ng
 xxran=mu;
 xold=xxran;
 ncmv=zeros(1,nr);
 xncmv=zeros(1,nr);
 xxncmv=zeros(1,nr);
 itte=0;
 uold=zeros(1,nr);
 maxitte=100;
 precisionbeta=1;
 lamda=0.5*ones(ng,1);
 coszetaold=ones(ng,1);
 sig=0;
 while (itte<maxitte&&precisionbeta>1e-6)
  itte=itte+1;
  xxold=xold;
  xxncmv=xncmv;
  xold=xxran;
  xncmv=ncmv;
  [g,gd,gx]=probconstraint(xval,xxran,i,ng,nr,ns);
  gnum=gnum+1;
  dgnum=dgnum+nr;
  gxx=gx(i,:);
  gdd=gd(i,1:n);
  [u,gu]=transform(xxran,gxx,type,mu,sigma); %转变成标准正态分布
  ncmv=-gu/norm(gu);
  kc=(ncmv-xncmv)*(xncmv-xxncmv).';
  if itte==1 
      coszeta=1;
  else
     coszeta=ncmv*xncmv'/norm(ncmv)/norm(xncmv);
  end
  if (sign(kc)<=0&&itte>2)  
       sig=1;
  end

  if sig==1 
      [lamdanew]=adpt(lamda(i),coszeta,coszetaold,betat);
      lamda(i)=lamdanew;
  end
  
  if sig==1   
      u=-betat*gu/norm(gu);
      udev=uold+lamda(i)*(u-uold);
      u=betat*udev/norm(udev);
  else
      u=-betat*gu/norm(gu);
  end
  uold=u;
  xxran=transforminv(u,type,mu,sigma);
  fvalx(i)=-g(i);
  precisionbeta=norm(xxran-xold)/norm(xxran);
 end
%  dgnum=dgnum+(n-ns);
 dfdx(i,:)=-gdd;
 xran(i,:)=xxran;
end

fval=fvalx';

%%%%%%%%%%%%%%%%%%%%%
%%%lamda
%%%%%%%%%%%%%%%%%%%%%
function [lamdanew]=adpt(lamdaold,coszeta,coszetaold,betat)
dlta=coszeta-coszetaold;
if dlta<0 %%%%
      ratiozeta=acos(coszetaold)/acos(coszeta); %%%%%角度的比值
      if ratiozeta>0.2
       lamdanew=lamdaold*acos(coszetaold)/acos(coszeta);
      else
       lamdanew=lamdaold*0.2;
      end
else
    lamdanew=lamdaold;
end
lamdanew=max([lamdanew,0.1]);


