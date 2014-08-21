function print_fsolve_out(fval, exitflag, output, jacobian, outFile)

[f_row,f_col] = size(fval);
fid1 = fopen(outFile,'a');

fprintf(fid1,'%s\n','fval:');

str_tmp = '\n';
for i = 1:f_col str_tmp =  ['%1.6f    ',str_tmp]; end;
for i = 1:f_row
    fprintf(fid1,str_tmp,fval(i,:));
end;
fprintf(fid1,'\n\n');
clear str_tmp;

fprintf(fid1,'%s','exitflag:'); fprintf(fid1,'%2d\n',exitflag);
if (exitflag<=0) fprintf(fid1,'%s\n\n','***********ERROR exitflag non positive!'); end;
%fprintf(fid1,'\n\n\n');

fprintf(fid1, '%s\n', 'output:');
fprintf(fid1, '%s   ','iterations:');       fprintf(fid1,'%1.0f\n',output.iterations); 
fprintf(fid1, '%s   ','funcCount:');        fprintf(fid1,'%1.0f\n',output.funcCount); 
fprintf(fid1, '%s   ','algorithm:');        fprintf(fid1,'%s\n',output.algorithm); 
fprintf(fid1, '%s   ','firstorderopt:');    fprintf(fid1,'%1.4e\n',output.firstorderopt); 
fprintf(fid1, '%s   ','message:');          fprintf(fid1,'%s\n\n\n',output.message); 

%{ 
printing jacobian
[jac_row,jac_col] = size(jacobian);
fprintf(fid1,'%s\n','jacobian:');

str_tmp = '\n';
for i = 1:jac_col str_tmp =  ['%1.4f    ',str_tmp]; end;
for i = 1:jac_row
    fprintf(fid1,str_tmp,jacobian(i,:));
end;
fprintf(fid1,'\n');
clear str_tmp;
%}

fprintf(fid1,'%s\n\n','----------------------------'); 
fclose(fid1);
end
