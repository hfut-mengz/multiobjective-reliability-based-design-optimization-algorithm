function [g,gd,gx]=probconstraint(xval,xran,i,ng,nr,ns)
g=fval(xval,xran);
n=length(xval); %%%��Ʊ����ĸ���

if n-ns>0    %%%%�����Ʊ������в�����������Ĳ���
 for i=1:1:n
     xnew=xval;
     if xnew(i)==0;
      xnew(i)=0.0001 ;
     else
      xnew(i)=xnew(i)*0.0001+xnew(i);
     end
     dltax=xnew(i)-xval(i);
     gnew=fval(xnew,xran);
     gd(:,i)=(gnew-g)/dltax;
 end
end
   nsi=n-ns+1;
   
 for i=1:1:nr
     xnew=xran;
     if xnew(i)==0
      xnew(i)=0.0001 ;
     else
      xnew(i)=xnew(i)*0.0001+xnew(i);
     end
     dltax=xnew(i)-xran(i);
     gnew=fval(xval,xnew);
     gx(:,i)=(gnew-g)/dltax;
 end
 
 nsi=n-ns+1;
 if ns>0
     gd(:,nsi:n)=gx(:,1:ns);
 end