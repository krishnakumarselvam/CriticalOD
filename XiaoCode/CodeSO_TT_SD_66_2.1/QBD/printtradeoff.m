function printtradeoff(tradeoff_py,tradeoff)
fid=fopen(tradeoff_py,'w');
fprintf(fid,'%s %f','tradeoff =',tradeoff);
fclose(fid);