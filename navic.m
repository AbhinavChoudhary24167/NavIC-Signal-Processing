for j = 1:4
    load(['in_data', num2str(j), '.mat']);  % Assumes variable name 'in'
    filename = sprintf('in_data%d.h', j);
    fid = fopen(filename, 'w');
    fprintf(fid, '#ifndef IN_DATA%d_H\n#define IN_DATA%d_H\n\n', j, j);
    fprintf(fid, 'double in_data%d[%d] = {\n', j, length(in));
    fprintf(fid, '%.10f,\n', in);
    fprintf(fid, '};\n\n#endif');
    fclose(fid);
end
