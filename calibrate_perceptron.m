function calibrate_perceptron(step_size, epoch_time, hidden_nodes_A_layer, hidden_nodes_B_layer, training_file, validating_file, testing_file, save_file_prefix)
%recommended step_size is 0.0003
%recommended epoch_time is 20000
%recommended hidden_nodes 3:5 for both layers
%training_file    : '/home/jessada/development/thesis/data/training/probabilistic_training_featuresset'
%validating_file  : '/home/jessada/development/thesis/data/training/probabilistic_validating_featuresset'
%testing_file     : '/home/jessada/development/thesis/data/training/probabilistic_testing_featuresset'
%save_file_prefix : '/home/jessada/development/thesis/data/training/params/training_probabilistic_2layer_params'

time_stamp                   = datestr(now, 'yyyymmddHHMMss');
result_dir                   = strcat(save_file_prefix, '/', time_stamp, '/');
mkdir(result_dir);

bias                         = [1 2];
sum_bias                     = sum(bias);

training_dataset             = importdata(training_file, '\t');
training_header              = training_dataset.rowheaders';
data_temp                    = training_dataset.data';
training_features            = data_temp(1:size(data_temp,1)-1, :);
training_targets             = data_temp(size(data_temp,1), :);
[unique_targets m n]         = unique(training_targets);
frequency                    = hist(training_targets(:),unique_targets) ./ bias;
training_targets_frequency   = reshape(frequency(n),size(training_targets));

validating_features_file     = importdata(validating_file, '\t');
validating_header            = validating_features_file.rowheaders';
data_temp                    = validating_features_file.data';
validating_features          = data_temp(1:size(data_temp,1)-1, :);
validating_targets           = data_temp(size(data_temp,1), :);
[unique_targets m n]         = unique(validating_targets);
frequency                    = hist(validating_targets(:),unique_targets) ./ bias;
validating_targets_frequency = reshape(frequency(n),size(validating_targets));

testing_features_file        = importdata(testing_file, '\t');
testing_header               = testing_features_file.rowheaders';
data_temp                    = testing_features_file.data';
testing_features             = data_temp(1:size(data_temp,1)-1, :);
testing_targets              = data_temp(size(data_temp,1), :);
[unique_targets m n]         = unique(testing_targets);
frequency                    = hist(testing_targets(:),unique_targets) ./ bias;
testing_targets_frequency    = reshape(frequency(n),size(testing_targets));

min_error          = 10;
model_error_matrix = [];
w1_struct_matrix   = {};
w2_struct_matrix   = {};
w3_struct_matrix   = {};
for i = hidden_nodes_A_layer
    model_error_array = [];
    w1_struct_array   = {};
    w2_struct_array   = {};
    w3_struct_array   = {};
    for j = hidden_nodes_B_layer
        disp(sprintf('run with hidden nodes A : %2d   B : %2d', i, j));
        drawnow;
        [weight_first_layer, weight_second_layer, weight_third_layer, error] = perceptron_learning(training_features, training_targets, training_targets_frequency, validating_features, validating_targets, validating_targets_frequency, step_size, epoch_time, i, j, result_dir, sum_bias); 
        model_error_array = [model_error_array, error];
        w1_struct_array   = [w1_struct_array, {weight_first_layer}];
        w2_struct_array   = [w2_struct_array, {weight_second_layer}];
        w3_struct_array   = [w3_struct_array, {weight_third_layer}];

        %Keep track of the minimum error rate for each number of hidden nodes
        if error < min_error
            min_error = error;
            best_weight_first_layer  = weight_first_layer;
            best_weight_second_layer = weight_second_layer;
            best_weight_third_layer  = weight_third_layer;
        end
    end;
    model_error_matrix = [model_error_matrix; model_error_array];
    w1_struct_matrix   = [w1_struct_matrix; w1_struct_array];
    w2_struct_matrix   = [w2_struct_matrix; w2_struct_array];
    w3_struct_matrix   = [w3_struct_matrix; w3_struct_array];
end;

%save weights from all configurations in case of possible future use
model_weight_matrix = struct('w1', w1_struct_matrix, 'w2', w2_struct_matrix, 'w3', w3_struct_matrix);

[out1, out2, out3] = forward_pass(best_weight_first_layer, best_weight_second_layer, best_weight_third_layer, testing_features);
biased_error = evaluate_model(out3, testing_targets, testing_targets_frequency, sum_bias);

%neutral_SNP_features = testing_features(:, find(testing_targets == 0));
%[neutral_out1, neutral_out2, neutral_out3] = forward_pass(best_weight_first_layer, best_weight_second_layer, best_weight_third_layer, neutral_SNP_features);

%pathogenic_SNP_features = testing_features(:, find(testing_targets == 1));
%[pathogenic_out1, pathogenic_out2, pathogenic_out3] = forward_pass(best_weight_first_layer, best_weight_second_layer, best_weight_third_layer, pathogenic_SNP_features);

f = figure;

mesh([hidden_nodes_B_layer], [hidden_nodes_A_layer], model_error_matrix);
xlabel('hidden nodes B');
ylabel('hidden nodes A');
zlabel('%Error');
title('Hidden nodes against Error rate');
model_error_figure_filename = strcat(result_dir, 'mesh.png');
print(f, '-dpng', model_error_figure_filename);

%{
neutral_distribution = hist(neutral_out3, [0:0.01:1]);
plot([0:0.01:1], neutral_distribution, 'b-');
hold on
pathogenic_distribution = hist(pathogenic_out3, [0:0.01:1]);
plot([0:0.01:1], pathogenic_distribution, 'r-');
legend('neutral probability distribution', 'pathogenic probability distribution');
hold off
ylabel('samples');
xlabel('pathogenic probability');
title(sprintf('probability distribution with biased error : %6.4f', biased_error));
probability_distribution_filename = strcat(result_dir, 'prob.png');
print(f, '-dpng', probability_distribution_filename);
%}
close(f);

params_file = sprintf('%s_%6.4f.mat', strcat(result_dir, '/combiner_params'), biased_error);
save(params_file);

plot_result(params_file, result_dir);
