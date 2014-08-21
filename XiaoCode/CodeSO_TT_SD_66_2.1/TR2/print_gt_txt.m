function print_gt_txt(fileTXT,fileData,fileNwk,fileQgParam,beforeRounding,gtNew)

load(fileData,'x_out');
load(fileNwk,'cycleTimes_s','phase2Node','availableGreenTimes_s');
load(fileQgParam,'last_index');

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

nbPhases = length(phase2Node);
if (beforeRounding==1)
    fprintf(fid,'GTs before rounding\n');
    currGT = x_out(last_index:(last_index+nbPhases-1)).*cycleTimes_s(phase2Node);
    if ~isempty(gtNew) error('before rounding the vector gtNew should be empty'); end;    
else
    fprintf(fid,'GTs after rounding\n');
    currGT = gtNew;
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
