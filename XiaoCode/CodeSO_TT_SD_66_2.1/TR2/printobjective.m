function printobjective(newAvgTTFile,SD)
% print the new objective function value SD*TT to the output file 
fid=fopen(newAvgTTFile,'r');

data = fscanf(fid,'%d %d %f %d');
fclose(fid);
newobj =  data(3)*SD;
newdata=[data(1),data(2),newobj,data(4)];
fid=fopen(newAvgTTFile,'w');
fprintf(fid,'%d %d %f %d',newdata);
fclose(fid);

