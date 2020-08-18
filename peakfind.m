function [pk,indpk,indv]=peakfind(h,k)
%h=input histogram(or homogram)
%pk: magnitude of peaks
%indpk: gray level assosiated with (the angle at which peak happened)
%indv: the angle at which valley happened
%k: if k=1 then the algorithm is compeletly implemented but if k=2 the
%phase of significant peak detection is ignored.

%finding all peaks:
d1=diff([h,0]);% 0 is added to consider last point
d2=[0.0001,d1];
d2(end)=[];
ind1=find((d2>0)&(d1<0));
pk1=h(ind1);
plot(ind1,pk1,'rx');
hold on;plot(h);

%significant peak detection finding peaks of peaks:
if k==1
    d1=diff([pk1,0]);% 0 is added to consider last point
    d2=[0.0001,d1];
    d2(end)=[];
    ind2=find((d2>0)&(d1<0));
    pk2=pk1(ind2);
    ind2=ind1(ind2);
else 
    ind2=ind1;
    pk2=pk1;
end


%step1: elimination of small peaks:
pkmax=max(pk2);
ind2(pk2/pkmax<0.02)=[];
pk2(pk2/pkmax<0.02)=[];
figure,plot(ind2,pk2,'go');

%step2: choose one peak if two peaks are too close
% d=diff(ind2);
% for i=1:length(d)
%     if d(i)<15
%         if pk2(i)>pk2(i+1)
%             ind2(i+1)=-10;
%             pk2(i+1)=-10;
%             try
%             d(i+1)=d(i+1)+d(i);
%             end
%         else
%             ind2(i)=-10;
%             pk2(i)=-10;
%         end
%     end
% end
% pk2(pk2==-10)=[];
% ind2(ind2==-10)=[];
% 
[pk2,kk]=sort(pk2,'descend');
ind2=ind2(kk);
for i=1:length(pk2)
    try
    pk2((ind2>ind2(i)-5)&(ind2<ind2(i)+5)&(ind2~=ind2(i)))=[];
    ind2((ind2>ind2(i)-5)&(ind2<ind2(i)+5)&(ind2~=ind2(i)))=[];
    catch
        break;
    end
end
[ind2,kk]=sort(ind2);
pk2=pk2(kk);
figure,plot(ind2,pk2,'k+');
    
%step3: remove peaks if the valley between peaks is not obvious:
ind=[];
for i=1:length(pk2)-1
    pav=(pk2(i)+pk2(i+1))/2;
    av=mean(h(ind2(i):ind2(i+1)));
    if (av>.75*pav)&(min(h(ind2(i):ind2(i+1)))>0.9*min(pk2(i),pk2(i+1)))
        if pk2(i)>pk2(i+1)
            ind=[ind,i+1];
        else
            ind=[ind,i];
        end
    end
end
pk2(ind)=[];
ind2(ind)=[];
% pk2((ind2<20)|(ind2>230))=[];
% ind2((ind2<20)|(ind2>230))=[];
%segmentation:
%II=255*ones(size(I))
for i=1:length(ind2)-1
    temp=ones(size(h))*1e10;
    temp(ind2(i):ind2(i+1))=h(ind2(i):ind2(i+1));
    [v(i),indv(i)]=min(temp);
    try
        m=indv(i-1);
    catch
        m=1;
    end
    %II((I>=m)&(I<indv(i)))=round(i*255/(length(ind2)));
end
indv; %thresholds
indpk=ind2;
pk=pk2;
plot(indpk,pk,'k')
%figure
%imshow(uint8(II))