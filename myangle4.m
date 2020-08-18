function [theta,dtheta,mrgim,label,indv,devind,labelcolim,labelim,corestab...
    ]=myangle4(I,bw,x,y,xnbr,ynbr,dist,trireg,tridevind,xtri,ytri,sig)
%% What does the function do? 
% This funxtion assign every pore an angle and draw a histogram and
% calculate the valleyes of histogram and cluster them and label the same ...
% directional angles(this function does not consider continuty  criteria and does not merge ...
% the same directional domains)
%% Input and output are as below:
% sig=0.001;
%modified histogram
%goals: assigning angels to pores, histogram analysis and clusering pores 
%based on the assigned angels
%J= labeled image;
%input params:
%   I=initail image
%   bw=bw image
%   x,y= pore coordinates
%   xnbr,ynbr=pore neighbors
%   dist= distance of each pore from its neighbors
%   trireg= a vector which indicates regular triangels by 1 and others by 0
%   tridevind= deviation of each triangel
%   xtri,ytri= vertices of each triangel.
%output params:
%   theta=assign angel between 0-60 to every center on the basis of its neighbors and assign 100 to disordered ones
%   dtheta=difference of theta with respect to its 6 neighbors that can be ...
%          used later for edge candidation (if we have few candidate beside each other)  and assign 100 to disordered ones
%   mrgim= final picture which is color coded classes merged with original image
%   devind= deviation for each pore(first we calculate the average of ...
%           dev_phi and dev_d for every triangle then we average these values between all 6 neighboring triangle 
%   label: every point is being assigned a label via angle direction criteria
%   indv: is a place of every valley e.g. 15 degree valley
%   tablecolim: it is a colored image which only has labels (no texture) and every label has a specific color
%   labelim: a picture which define labels with numbers and when we show it ...
%            it is a black and white picture which shows labeled area through white color
%   corestab: relates index of triangle to the index of x's and define ...
%             which points area related to which triangle and it is in contrast with ind
%__________________________________________________________________________
%% determining the correspondance between triangle list and pore list
corestab=zeros(6,length(x));%determines the neighboring triangel inidices for a pore.
for i=1:length(x)
    indx=find(xtri(1,:)==x(i));%x is a center of a pore and (xtri(1,:) is a one vertex of a triangle
    indy=find(ytri(1,:)==y(i));
    indxy1=intersect(indx,indy);% define common triangles
    indx=find(xtri(2,:)==x(i));
    indy=find(ytri(2,:)==y(i));
    indxy2=intersect(indx,indy);
    indx=find(xtri(3,:)==x(i));
    indy=find(ytri(3,:)==y(i));
    indxy3=intersect(indx,indy);
    %% Fantastic
    indxy=union(union(indxy1,indxy2),indxy3);
    corestab(1:length(indxy),i)=indxy;
end
%__________________________________________________________________________
%% assigning angels:
[n1,n2]=size(bw);
dtheta=100*ones(6,length(x)); % first suppose every triangle is disordered
devind=100*ones(size(x));
for i=1:length(x)
    try
    d=dist(:,i);
    xn=xnbr(:,i);
    yn=ynbr(:,i);
    d(d==0)=[]; %null the invalid distances
    xn(xn==0)=[];
    yn(yn==0)=[];
    xd=xn-x(i);
    yd=yn-y(i);
    phi=atan2(yd,xd)*180/pi;
    phi1=sort(phi);
    phi1(end+1)=phi1(1)+360;
    %if (length(phi)==6)&&(min(diff(phi1))>30)&&(max(diff(phi1))<85)
    %check for regularity of neighboring triangels:
    cond=(1==1);
    devind(i)=0;
    tmp=corestab(:,i);
    tmp(tmp==0)=[];% e.g. we may have 10 triangle for a pore
    if length(tmp)==6 %if choose ordered areas (first criteria for ordered area)
    for j=1:size(corestab,1)
        if corestab(j,i)~=0
            cond=cond&&trireg(corestab(j,i));
            devind(i)=devind(i)+1/6*tridevind(corestab(j,i));%devind is for every point
        end
    end
    else 
        devind(i)=100;%greatly distorted.
    end
    
    if (cond)&&(length(xn)==6)
        dtheta(:,i)=diff(phi1);
        t=mod(phi,60);
        t1=t(t<30); % we want not to have misinterpretation about the 30 and 60(or 0)
        t2=t(t>=30);
        if mean(t2)-mean(t1)>30
            th1=mean([t1;t2-60]);
            th2=mean([t1+60;t2]);
            if th1>=0 %result must be between 0 and 59
                theta(i)=th1;
            else
                theta(i)=th2;
            end
        else
            theta(i)=mean(t); %normal state and normal averaging
        end
    else
        theta(i)=100;
        devind(i)=100;%greatly distorted
    end
    end
end
%__________________________________________________________________________
%% histogram analysis and clustering:
theta1=theta;
theta1(theta>=60-0.0001)=60-0.0001;%????????????????????
% theta1(theta==100)=[];
% h=hist(theta1,60);
%% weighted histogram
den=0;
hm=zeros(1,60);

for i=1:length(x)

    if theta1(i)~=100
%         sig=0.0010;
        hm(floor(theta1(i))+1)=hm(floor(theta1(i))+1)+exp(-(devind(i)).^2/sig);% sig is empirical
        den=den+exp(-(devind(i)).^2/sig);
    end

end
hm=hm/den;

e=(conv(hm,hamming(5))/sum(hamming(5))); %Smooth the histogram with a hamming window
e(1:2)=[];e(end-1:end)=[]; %Hamming has two extra point which must be eliminated
figure,plot(e,'r');
e=[e,e]; %attaching 60 to 0 degree
%% Peak fine function
[pk,indpk,indv]=peakfind(e,2);%modify peakfind circularly
%pk: magnitude of peaks
%indpk: the angle at which peak happened
%indv: the angle at which valley happened
%% Labeling
indv(indv>60)=[];
indv=[indv,61,inf];%????????????
% J=zeros(size(bw));
label=zeros(size(x));
for i=1:length(theta)%length(theta) is equal to length(x)
%     J(max(x(i)-4,1),max(y(i)-3,1):min(y(i)+3,n2))=min(find(theta(i)<indv));
%     J(max(x(i)-3,1),max(y(i)-4,1):min(y(i)+4,n2))=min(find(theta(i)<indv));
%     J(max(x(i)-2,1):min(x(i)+2,n1),max(y(i)-5,1):min(y(i)+5,n2))=min(find(theta(i)<indv));
%     J(min(x(i)+4,n1),max(y(i)-3,1):min(y(i)+3,n2))=min(find(theta(i)<indv));
    label(i)=min(find(theta(i)<indv)); %assign clusters to the angles
end

% J(J==length(indv)-1)=1;
label(label==length(indv)-1)=1; %the last and the first label may be equall if the angle is not 60 (61) so the last and the first are the same
% J(J==length(indv))=length(indv)-1;
label(label==length(indv))=0; %the disordered area has a 0 labelwhich has angel between 61 and inf (their value is 100)
%% Attaching the label with texture
[mrgim,labelcolim,labelim]=imerg(I,label,x,y,xnbr,ynbr);
% figure;
% % imshow(label2rgb(J));
% J(J==numel(indv)-1)=0;
% im1=label2rgb(J);
% im1=rgb2hsv(im1);
% if size(I,3)==3
%     I=rgb2gray(I);
% end
% im1(:,:,3)=double(I)/255.*im1(:,:,3);
% im1=hsv2rgb(im1);
% imshow(im1)
% 
% % I(:,:,2)=40*J;
% % I(:,:,1)=30*(5-J)
% % imshow(I)
%% Execution of the program
% imshow(uint8(labelcolim))
% imshow((labelim))
% imshow(mrgim)

