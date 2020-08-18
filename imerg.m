function [mrgim,labelcolim,labelim]=imerg(I,label,x,y,xnbr,ynbr);
%% What does the function do? 
%merge the label with texture
%% Programs
if size(I,3)==3
    I=rgb2gray(I);
end
labelim=zeros(size(I));
for i=1:length(label)%the type of label e.g. is 5 but the length(label) is equal to length(x)
    for n1=x(i)-5:x(i)+5 %5 is changed to 7
        for n2=y(i)-5:y(i)+5
            try
            d1=sqrt((x(i)-n1)^2+(y(i)-n2)^2);%d1 has 1parameter
            d2=sqrt((xnbr(:,i)-n1).^2+(ynbr(:,i)-n2).^2);%d2 has 20 parameter
            cnd=(d2>=d1);
            cnd(xnbr(:,i)==0)=[];%some of the 20 are invalid
            I(n1,n2);
            if all(cnd)%if the point is inside the hexagonalit givesthe color of that hexagonal otherwise nothing happened
                labelim(n1,n2)=label(i);
            end
            end
        end
    end
end


%% Fantastic power 2
%% For changing the color very smoothly and periodic changing (e.g. from 0 to 59)
m=200;
Kr=m*trimf(labelim/max(max(label)),[0,.25,1])+256-m;
Kg=m*trimf(labelim/max(max(label)),[0,.75,1])+256-m;
Kb=m*max(trimf(labelim/max(max(label)),[-1,0,.5]),...
    trimf(labelim/max(max(label)),[.5,1,2]))+256-m;% for two part trimf
Kr(labelim==0)=255;
Kg(labelim==0)=255;
Kb(labelim==0)=255;%disordered area is whited
labelcolim(:,:,1)=Kr;labelcolim(:,:,2)=Kg;labelcolim(:,:,3)=Kb;
I4(:,:,1)=I;I4(:,:,2)=I;I4(:,:,3)=I;
I4=double(I4);
mrgim=uint8(I4/255.*labelcolim);
figure;imshow(mrgim);

