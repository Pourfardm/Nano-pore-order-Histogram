function [nlabel,mrgim,labelcolim,labelim]=myspreading(I,label,x,y,xnbr,ynbr);
%Connect neighboring domain to each other and give the the same label
%label: old label
%nlabel: new label
%% spreading algorithm based on our approach
%I=I1;
clabel=zeros(size(label));%class label

lab=1; %lab show the number of new label
for j=1:max(label)
    prog=not(label==j); %we only consider label j not the other label
    %indicates to the investigation progress, 
    %investigated triangels are labeled one.
ind=find(prog==0);%indicator of triangles to be considered
ind=ind(1);
nlabel(ind(1))=lab;
prog(ind(1))=1;
while ~isempty(ind)
%% Fantastic 3  
%     indx=find(abs(xc-xc(ind(1)))<=15);
%     indy=find(abs(yc-yc(ind(1)))<=15);
%     indxy=intersect(indx,indy);
    for i=1:6;
        indx=find(x==xnbr(i,ind(1)));
        indy=find(y==ynbr(i,ind(1)));
        indxy=intersect(indx,indy);%calculates the coordinates of neighbors of a pixel
        if (~isempty(indxy))&&(label(indxy)==label(ind(1)))&&(prog(indxy)==0);%if the label of neighbor equall to center pixel and if it does not investigated before
            ind(end+1)=indxy;
            nlabel(indxy)=lab;
            prog(indxy)=1;
        end
    end
    ind(1)=[]; %now we just only consider one center
    if isempty(ind)
        ind=find(prog==0);%indicator of triangles to be considered
        if ~isempty(ind)
            ind=ind(1);
            lab=lab+1
            nlabel(ind(1))=lab;
            prog(ind(1))=1;
%             I(round(xc(indxy(i))),round(yc(indxy(i))),1)=255;
%             I(round(xc(indxy(i))),round(yc(indxy(i))),2)=0;
%             I(round(xc(indxy(i))),round(yc(indxy(i))),3)=0;
%             imshow(I)
%             pause
        end
        h=hist(prog);
        h(1);
    end
end
end
%elimination of one triangle regions:
% lab=randperm(lab);
% nnlabel=nlabel;
%  for i=1:length(lab)
% %      if sum(trilabel==i)==1
%          ntrilabel(trilabel==i)=lab(i);
% %      end
%  end
% ind=sort(unique(trilabel));
% for i=0:length(ind)-1
%     trilabel(trilabel==ind(i+1))=i;
% end
    
[mrgim,labelcolim,labelim]=imerg(I,label,x,y,xnbr,ynbr);        
[mrgim,labelcolim,labelim]=imerg(I,nlabel,x,y,xnbr,ynbr);%random color would be better?????????????????????????
