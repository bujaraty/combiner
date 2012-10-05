function plot_result(params_file, output_dir)
%recent param 'params_file' : '/home/jessada/development/thesis/data/training/params/selected/training_probabilistic_params_0.3530_20120326162053.mat'

load(params_file);

neutral_SNP_index    = find(test_targets == 0);
neutral_SNP_features = test_features(:, neutral_SNP_index);
[number_of_features number_of_neutral_samples] = size(neutral_SNP_features);

[neutral_out1, neutral_out2, neutral_out3] = forward_pass(best_weight_first_layer, best_weight_second_layer, best_weight_third_layer, neutral_SNP_features);

pathogenic_SNP_index    = find(test_targets == 1);
pathogenic_SNP_features = test_features(:, pathogenic_SNP_index);
number_of_pathogenic_samples = size(pathogenic_SNP_features, 2);

[pathogenic_out1, pathogenic_out2, pathogenic_out3] = forward_pass(best_weight_first_layer, best_weight_second_layer, best_weight_third_layer, pathogenic_SNP_features);

min_gerp = min(neutral_SNP_features(1,:));
max_gerp = max(neutral_SNP_features(1,:));

gerp_roc_range = min_gerp:(max_gerp-min_gerp)/1000:max_gerp;
FP_rate =[];
TP_rate = [];

for gerp_pathogenic_ratio = gerp_roc_range
    FP_rate = [FP_rate, sum([neutral_SNP_features(1,:)] > gerp_pathogenic_ratio, 2) / number_of_neutral_samples];
    TP_rate = [TP_rate, sum([pathogenic_SNP_features(1,:)] > gerp_pathogenic_ratio, 2) / number_of_pathogenic_samples];
end;

roc_range = 0.0:0.001:1;
FP_rates = [];
TP_rates = [];

for pathogenic_ratio = roc_range
    FP_rates = [FP_rates, sum([neutral_SNP_features(2:number_of_features,:); neutral_out3] > pathogenic_ratio, 2) / number_of_neutral_samples];
    TP_rates = [TP_rates, sum([pathogenic_SNP_features(2:number_of_features,:); pathogenic_out3] > pathogenic_ratio, 2) / number_of_pathogenic_samples];
end;

FP_rates = [FP_rate; FP_rates];
TP_rates = [TP_rate; TP_rates];

f = figure;
hold on
plot(FP_rates(1,:), TP_rates(1,:), '--m');
plot(FP_rates(2,:), TP_rates(2,:), '--b');
plot(FP_rates(3,:), TP_rates(3,:), '--g');
plot(FP_rates(4,:), TP_rates(4,:), '--k');
plot(FP_rates(5,:), TP_rates(5,:), '--c');
plot(FP_rates(6,:), TP_rates(6,:), '--r');
plot(FP_rates(7,:), TP_rates(7,:), 'k');
hold off
ylabel('True positive rate');
xlabel('False positive rate');
title('ROC curve');
legend('GERP++', 'SIFT', 'polyphen2', 'phylop', 'LRT', 'Mutation Taster', 'Combiner');
roc_filename = strcat(output_dir, 'roc.png');
print(f, '-dpng', roc_filename);
%{
subplot(2,3,1);
neutral_distribution = hist(neutral_SNP_features(2,:), [0:0.02:1]);
plot([0:0.02:1], neutral_distribution, '-.b');
hold on
pathogenic_distribution = hist(pathogenic_SNP_features(2,:), [0:0.02:1]);
plot([0:0.02:1], pathogenic_distribution, '-.r');
legend('neutral variants score', 'pathogenic variants score');
hold off
ylabel('samples');
xlabel('score');
title('SIFT score distribution');
%}
subplot(2,2,1);
neutral_distribution = hist(neutral_SNP_features(3,:), [0:0.02:1]);
plot([0:0.02:1], neutral_distribution, '-.b');
hold on
pathogenic_distribution = hist(pathogenic_SNP_features(3,:), [0:0.02:1]);
plot([0:0.02:1], pathogenic_distribution, '-.r');
legend('neutral variants score', 'pathogenic variants score');
hold off
ylabel('samples');
xlabel('score');
title('Polyphen2 score distribution');
%{
subplot(2,3,3);
neutral_distribution = hist(neutral_SNP_features(4,:), [0:0.02:1]);
plot([0:0.02:1], neutral_distribution, '-.b');
hold on
pathogenic_distribution = hist(pathogenic_SNP_features(4,:), [0:0.02:1]);
plot([0:0.02:1], pathogenic_distribution, '-.r');
legend('neutral variants score', 'pathogenic variants score');
hold off
ylabel('samples');
xlabel('score');
title('Phylop score distribution');
%}
subplot(2,2,2);
neutral_distribution = hist(neutral_SNP_features(6,:), [0:0.02:1]);
plot([0:0.02:1], neutral_distribution, '-.b');
hold on
pathogenic_distribution = hist(pathogenic_SNP_features(6,:), [0:0.02:1]);
plot([0:0.02:1], pathogenic_distribution, '-.r');
legend('neutral variants score', 'pathogenic variants score');
hold off
ylabel('samples');
xlabel('score');
title('Mutation Taster score distribution');

subplot(2,2,3);
neutral_distribution = hist(neutral_SNP_features(1,:), min_gerp:(max_gerp-min_gerp)/50:max_gerp);
plot([0:0.02:1], neutral_distribution, '-.b');
hold on
pathogenic_distribution = hist(pathogenic_SNP_features(1,:), min_gerp:(max_gerp-min_gerp)/50:max_gerp);
plot([0:0.02:1], pathogenic_distribution, '-.r');
legend('neutral variants score', 'pathogenic variants score');
hold off
ylabel('samples');
xlabel('score');
title('GERP++ score distribution');

subplot(2,2,4);
neutral_distribution = hist(neutral_out3, [0:0.02:1]);
plot([0:0.02:1], neutral_distribution, '-.b');
hold on
pathogenic_distribution = hist(pathogenic_out3, [0:0.02:1]);
plot([0:0.02:1], pathogenic_distribution, '-.r');
legend('neutral variants score', 'pathogenic variants score');
hold off
ylabel('samples');
xlabel('score');
title('Combiner score distribution');
score_distribution_filename = strcat(output_dir, 'score_distribution_combiner.png');
print(f, '-dpng', score_distribution_filename);

close(f);
