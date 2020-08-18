function [xtri,ytri,xc,yc,trireg,tridevind]=triangulation_ab(I,x,y,xnbr,ynbr,tol_phi,tol_d)
%goals: finding triangels and order determination, 

%inputs: x,y,xnbr,ynbr,
%finding xc,yc cordinates of triangle centroids;
%finding xrti,ytri cordinates of triangle vertices for all triangels.
%test for regularity => regular ones are indexd 1 in trireg
%tridevind= regularity index for each of triangels=sum of deviations.
%--------------------------------------------------------------------------
%finding triangels:
n=length(x);
xtri=[];
ytri=[];
for i=1:n
    try
    xn=xnbr(:,i);
    yn=ynbr(:,i);
    xn(xn==0)=[];
    yn(yn==0)=[];
    phi=atan2((xn-x(i)),(yn-y(i)))*180/pi;
    [phi,k]=sort(phi);
    xn=xn(k);xn(end+1)=xn(1);
    yn=yn(k);yn(end+1)=yn(1);
    for j=1:length(xn)-1
        [xtri(:,end+1),k]=sort([x(i),xn(j),xn(j+1)]);
        yk=[y(i),yn(j),yn(j+1)];
        ytri(:,end+1)=yk(k);
    end
    catch
        1
    end
end
xc=mean(xtri);
yc=mean(ytri);
[a,ind]=unique([xc;yc]','rows');
a=a';
xc=a(1,:);
yc=a(2,:);
xtri=xtri(:,ind);
ytri=ytri(:,ind);
%--------------------------------------------------------------------------
%triangular regularity test:
phi=zeros(1,3);
trireg=zeros(size(xc));
% tol_phi=.15;
% tol_d=.15;
for i=1:length(xc)
    d1=sqrt((xtri(1,i)-xtri(2,i))^2+(ytri(1,i)-ytri(2,i))^2);
    d2=sqrt((xtri(1,i)-xtri(3,i))^2+(ytri(1,i)-ytri(3,i))^2);
    d3=sqrt((xtri(3,i)-xtri(2,i))^2+(ytri(3,i)-ytri(2,i))^2);
    dm=(d1+d2+d2)/3;
    dev_d=sqrt(((d1-dm)^2+(d2-dm)^2+(d3-dm)^2)/3)/dm;
    v12=[xtri(2,i)-xtri(1,i),ytri(2,i)-ytri(1,i)];
    v12=v12/sqrt((v12(1))^2+(v12(2))^2);
    v23=[xtri(3,i)-xtri(2,i),ytri(3,i)-ytri(2,i)];
    v23=v23/sqrt((v23(1))^2+(v23(2))^2);
    v31=[xtri(1,i)-xtri(3,i),ytri(1,i)-ytri(3,i)];
    v31=v31/sqrt((v31(1))^2+(v31(2))^2);
    phi(1)=acos(v12*(-v31)')*180/pi;
    phi(2)=acos(v23*(-v12)')*180/pi;
    phi(3)=acos(v31*(-v23)')*180/pi;
    dev_phi=sqrt(((phi(1)-60)^2+(phi(2)-60)^2+...
        (phi(3)-60)^2)/3)/60;
    trireg(i)=(dev_phi<=tol_phi)&&(dev_d<=tol_d);
    tridevind(i)=0.5*(dev_phi+dev_d);
end


%[mrgim,labelcolim]=triimerg(I,trireg,xtri,ytri);