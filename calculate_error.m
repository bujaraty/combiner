function [error] = calculate_error(output, targets, targets_frequency, sum_bias)
error = (sum(sum(abs(output - targets) ./ targets_frequency))) / sum_bias;
