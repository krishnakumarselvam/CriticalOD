function print_gt_txt_notFromQg(fileTXT,fileNwk,beforeRounding,gtVect)
% gtVect :  if beforeRounding==1 :   gtInit is the equivalent of x_out(last_index:(last_index+nbPhases-1))
%                                   it will be multiplied by cycleTimes_s(phase2Node) to be converted to actual green times
%           else it will be the actual rounded green times  (no multiplication by cycle times)


load(fileNwk,'cycleTimes_s','phase2Node','availableGreenTimes_s');
nbPhases = length(phase2Node);

%file_p='';

if ~strcmp(fileTXT,'')
    if (beforeRounding==1)
        fid=fopen(fileTXT,'w');
    else
        fid=fopen(fileTXT,'a');
    end;
else fid=1;
end;

fprintf(fid,'\n************\n\n');

if (beforeRounding==1)
    fprintf(fid,'GTs before rounding\n');
    currGT = gtVect.*cycleTimes_s(phase2Node);
else
    fprintf(fid,'GTs after rounding\n');
    currGT = gtVect;
end;    

fprintf(fid,'NodeIndx availableGT\n');
for p = 1:length(unique(phase2Node))
    fprintf(fid,'%d %3.4f\n',p,availableGreenTimes_s(p));
end;

fprintf(fid,'\n\n');

fprintf(fid,'NodeIndx phaseGT\n');
for p = 1:nbPhases
    fprintf(fid,'%d %2.4f\n',[phase2Node(p),currGT(p)]);
end;

if ~strcmp(fileTXT,'')
    fclose(fid);
end;
