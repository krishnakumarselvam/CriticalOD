function SD = refine(entrTTfile,exitTTfile);


% C142 = [folderName,'Centroid142.txt'];
% C143 = [folderName,'Centroid143.txt'];
% C759 = [folderName,'Centroid759.txt'];
% C760 = [folderName,'Centroid760.txt'];


fpp=fopen(exitTTfile,'r');
tem2=fscanf(fpp,'%f');
fclose(fpp);
tem2=reshape(tem2,3,length(tem2)/3)';
exitT=unique(tem2,'rows');
exitTasc = sortrows(exitT,[1 2]);

fp=fopen(entrTTfile,'r');
tem1=fscanf(fp,'%f');
fclose(fp);
tem1=reshape(tem1,3,length(tem1)/3)';
entrT=unique(tem1,'rows');
entrTasc = sortrows(entrT,[1 2]);

%for those vehicles entered from 759£¬760,759¡¢760 is both entrance and exit;travel
%time for them will not be considered.
% fp=fopen(C759,'r');%%those vehicle entered from 759
% ID759=fscanf(fp,'%d');
% fclose(fp);
% q=length(ID759);
% A=find(exitTasc(:,3)==759);
% for i=1:q;
%     for j=1:length(A);
%         if exitTasc(A(j),1)==ID759(i);
%             exitTascdel759(i)=A(j);
%         end;
%     end;
% end;
% exitTasc(exitTascdel759,:)=[];%%delete those vehicles entered and exit from 759
% 
% 
% fp=fopen(C760,'r');
% ID760=fscanf(fp,'%d');
% fclose(fp);
% q=length(ID760);
% A=find(exitTasc(:,3)==760);
% for i=1:q;
%     for j=1:length(A);
%         if exitTasc(A(j),1)==ID760(i);
%             exitTascdel760(i)=A(j);
%         end;
%     end;
% end;
% exitTasc(exitTascdel760,:)=[];
% 
% fp=fopen(C143,'r');
% ID143=fscanf(fp,'%d');
% fclose(fp);
% q=length(ID143);
% A=find(exitTasc(:,3)==143);
% for i=1:q;
%     for j=1:length(A);
%         if exitTasc(A(j),1)==ID143(i);
%             exitTascundel143(i)=A(j);
%         end;
%     end;
% end;
% B=setdiff(A,exitTascundel143');
% exitTasc(B,:)=[];
% 
% fp=fopen(C142,'r');
% ID142=fscanf(fp,'%d');
% fclose(fp);
% q=length(ID142);
% A=find(entrTasc(:,3)==142);
% for i=1:q;
%     for j=1:length(A);
%         if entrTasc(A(j),1)==ID142(i);
%             entrTascundel142(i)=A(j);
%         end;
%     end;
% end;
% C=setdiff(A,entrTascundel142');
% entrTasc(C,:)=[];

IterExit=length(exitTasc);
for i = 1:IterExit;
    Iterexit{exitTasc(i,1)}=find(exitTasc(:,1)==exitTasc(i,1));
    lengthIterexit{exitTasc(i,1)}=length(Iterexit{exitTasc(i,1)});
    exitall(exitTasc(i,1),1:length(Iterexit{exitTasc(i,1)}))=Iterexit{exitTasc(i,1)};
    exitallTT(exitTasc(i,1),1:length(Iterexit{exitTasc(i,1)}))=exitTasc(exitall(exitTasc(i,1),1:length(Iterexit{exitTasc(i,1)})),2);
end

IterEntr=length(entrTasc);
for i = 1:IterEntr;
    Iterentr{entrTasc(i,1)}=find(entrTasc(:,1)==entrTasc(i,1));
    
    lengthIterentr{entrTasc(i,1)}=length(Iterentr{entrTasc(i,1)});
    entrall(entrTasc(i,1),1:length(Iterentr{entrTasc(i,1)}))=Iterentr{entrTasc(i,1)};
    entrallTT(entrTasc(i,1),1:length(Iterentr{entrTasc(i,1)}))=entrTasc(entrall(entrTasc(i,1),1:length(Iterentr{entrTasc(i,1)})),2);
end

[rexit,cexit]=size(exitall);
[rentr,centr]=size(entrall);

if centr>cexit;
    entrallTT(:,(cexit+1):centr)=[];
else
    if centr<cexit;
        error(' nb of exited vehicles cannot exceed nb of entered vehicles');
    end;
end;

if rentr>rexit;
    entrallTT((rexit+1):rentr,:)=[];
else
    if rentr<rexit;
        error(' nb of exited vehicles cannot exceed nb of entered vehicles');
    end;
end;
% [REXIT,CEXIT]=find(exitallTT>0);
% AA=entrallTT(REXIT,CEXIT);
% [R,C]=find(AA==0);

travelT=exitallTT-entrallTT;
%calculate expected travel time
Z=find(travelT<0);
travelT(Z)=0;
nbvehicle=sum(travelT(:)~=0);
TT=sum(travelT(:));
EE=TT/nbvehicle/60;


TTID=find(travelT(:)>0);
varvehicle = var(travelT(TTID));
SD=varvehicle^0.5/60;

delete(entrTTfile);
delete(exitTTfile);




