function [best_weight_first_layer, best_weight_second_layer, best_weight_third_layer, min_error] = perceptron_learning(training_data, training_targets, training_targets_frequency, validating_data, validating_targets, validating_targets_frequency, step_size, epoch_time, hidden_nodes_A_layer, hidden_nodes_B_layer, figure_file_prefix, sum_bias)

alpha = 0.9;

[feature_size, training_datasize] = size(training_data);

min_error = 10;

weight_first_layer  = 0.01 .* randn(hidden_nodes_A_layer, feature_size + 1);
weight_second_layer = 0.01 .* randn(hidden_nodes_B_layer, hidden_nodes_A_layer + 1);
weight_third_layer  = 0.01 .* randn(1, hidden_nodes_B_layer + 1);

delta_weight_first_layer  = zeros(hidden_nodes_A_layer, feature_size + 1);
delta_weight_second_layer = zeros(hidden_nodes_B_layer, hidden_nodes_A_layer + 1);
delta_weight_third_layer  = zeros(1, hidden_nodes_B_layer + 1);

backprop_error   = [];
validating_error = [];

for i = 1:epoch_time
    %Forward Pass
    [out1, out2, out3] = forward_pass(weight_first_layer, weight_second_layer, weight_third_layer, training_data);
    
    %Backward Pass
    %size((out3 - training_targets))
    %a = (out3 - training_targets);
    %a(1)
    %b = (out3 - training_targets)./training_targets_frequency;
    %b(1)
    delta_o = (((out3 - training_targets) * min(training_targets_frequency)) ./ training_targets_frequency) .* ((1 + out3) .* (1 - out3)) * 0.25;
    delta_h = (weight_third_layer' * delta_o) .* ((1 + out2) .* (1 - out2)) * 0.5;
    delta_h = delta_h(1:hidden_nodes_B_layer, :);    
    delta_g = (weight_second_layer' * delta_h) .* ((1 + out1) .* (1 - out1)) * 0.5;
    delta_g = delta_g(1:hidden_nodes_A_layer, :);    
    
    delta_weight_first_layer  = (delta_weight_first_layer .* alpha) - (delta_g * [training_data ; ones(1, training_datasize)]') .* (1 - alpha);
    delta_weight_second_layer = (delta_weight_second_layer .* alpha) - (delta_h * out1') .* (1-alpha);
    delta_weight_third_layer  = (delta_weight_third_layer .* alpha) - (delta_o * out2') .* (1-alpha);   
    weight_first_layer  = weight_first_layer + delta_weight_first_layer .* step_size;
    weight_second_layer = weight_second_layer + delta_weight_second_layer .* step_size;
    weight_third_layer  = weight_third_layer + delta_weight_third_layer .* step_size;

    error              = calculate_error(out3, training_targets, training_targets_frequency, sum_bias);
    backprop_error     = [backprop_error, error];
    
    %Validating
    [out1, out2, out3] = forward_pass(weight_first_layer, weight_second_layer, weight_third_layer, validating_data);

    error              = calculate_error(out3, validating_targets, validating_targets_frequency, sum_bias);
    validating_error   = [validating_error, error];
    
    %Keep track of minimum error rate
    if error < min_error
        min_error = error;
        best_weight_first_layer  = weight_first_layer;
        best_weight_second_layer = weight_second_layer;
        best_weight_third_layer  = weight_third_layer;
    end
end;

f = figure;
plot(validating_error, 'b-');
hold on
plot(backprop_error, 'r-');
legend('Validation Error', 'BackProp Error');
hold off
ylabel('%Error');
xlabel('Epoch');
filename = sprintf('%splot_%02d_%02d.png', figure_file_prefix, hidden_nodes_A_layer, hidden_nodes_B_layer);
title('Number of iteration VS Error rate');
print(f, '-dpng', filename);
close(f);


