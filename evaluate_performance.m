function predict_pathogenic_probability(params_file, output_file)

load(params_file);


[out1, out2, out3] = forward_pass(best_weight_first_layer, best_weight_second_layer, best_weight_third_layer, test_features);
[pathstr,name,ext] = fileparts(output_file);

content_file = strcat(pathstr, '/content.txt');
header_file  = strcat(pathstr, '/header.txt');
dlmwrite(content_file, [round(out3*10000)/10000]', '\t');
dlmwrite(header_file, [char(test_header')], '');

tmp_file  = strcat(pathstr, '/tmp.txt');
command = sprintf('paste %s %s | sort -k1 > %s', header_file, content_file, tmp_file);
system([command]);
command = sprintf('sed ''s/ \\t/\\t/g'' %s > %s', tmp_file, output_file);
system([command]);
command = sprintf('rm %s ', header_file);
system([command]);
command = sprintf('rm %s ', content_file);
system([command]);
command = sprintf('rm %s ', tmp_file);
system([command]);

disp(sprintf('Predicted result is at : %s', output_file));


