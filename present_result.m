function present_result(params_file)
%recent param 'params_file' : '/home/jessada/development/thesis/data/training/params/selected/training_probabilistic_params_0.3530_20120326162053.mat'
%recent param 'testing_file'   : '/home/jessada/development/thesis/data/training/probabilistic_testing_dataset_3646'

load(params_file);

pathogenic_ratio     = 0.5;

%{
tmp_file             = importdata(testing_file,'\t');
data_temp            = tmp_file.data';
testing_info            = tmp_file.rowheaders';
testing_features        = data_temp(1:size(data_temp,1)-1, :);
testing_targets         = data_temp(size(data_temp,1), :);
%}
neutral_SNP_index    = find(testing_targets == 0);

%neutral_SNP_info     = testing_info(:, neutral_SNP_index);
neutral_SNP_features = testing_features(:, neutral_SNP_index);
neutral_SNP_targets  = testing_targets(:, neutral_SNP_index);


[out1, out2, out3] = forward_pass(best_weight_first_layer, best_weight_second_layer, best_weight_third_layer, neutral_SNP_features);
%neutral_error            = evaluate_probabilistic_model(out2, neutral_SNP_targets, size(neutral_SNP_info, 2) .* ones(size(neutral_SNP_targets)));
neutral_error = (sum(sum(abs(out3 - neutral_SNP_targets) ./ (size(neutral_SNP_features, 2) .* ones(size(neutral_SNP_targets))))));

false_positive_group = find(out3 > pathogenic_ratio);
true_negative_group  = find(out3 <= pathogenic_ratio);

number_of_neutral_SNP    = size(neutral_SNP_features, 2);

number_of_false_positive = size(false_positive_group, 2);
number_of_true_negative  = size(true_negative_group, 2);

%------------------------------- reserved ---------------------------------
%{
false_positive_neutral_SNP_info     = neutral_SNP_info(false_positive_group);
false_positive_neutral_SNP_features = neutral_SNP_features(:, false_positive_group);
false_positive_neutral_SNP_targets  = neutral_SNP_targets(:, false_positive_group);

true_negative_neutral_SNP_info      = neutral_SNP_info(true_negative_group);
true_negative_neutral_SNP_features  = neutral_SNP_features(:, true_negative_group);
true_negative_neutral_SNP_targets   = neutral_SNP_targets(:, true_negative_group);

[out1, out2] = forward_probabilistic_pass_2layer(best_weight_first_layer, best_weight_second_layer, false_positive_neutral_SNP_features);
error = evaluate_probabilistic_model(out2, false_positive_neutral_SNP_targets, size(false_positive_neutral_SNP_info, 2) .* ones(size(false_positive_neutral_SNP_targets)))

[out1, out2] = forward_probabilistic_pass_2layer(best_weight_first_layer, best_weight_second_layer, true_negative_neutral_SNP_features);
error = evaluate_probabilistic_model(out2, true_negative_neutral_SNP_targets, size(true_negative_neutral_SNP_info, 2) .* ones(size(true_negative_neutral_SNP_targets)))
%}
%------------------------------- reserved ---------------------------------

pathogenic_SNP_index    = find(testing_targets == 1);

%pathogenic_SNP_info     = testing_info(:, pathogenic_SNP_index);
pathogenic_SNP_features = testing_features(:, pathogenic_SNP_index);
pathogenic_SNP_targets  = testing_targets(:, pathogenic_SNP_index);

[out1, out2, out3] = forward_pass(best_weight_first_layer, best_weight_second_layer, best_weight_third_layer, pathogenic_SNP_features);
%pathogenic_error            = evaluate_probabilistic_model(out2, pathogenic_SNP_targets, size(pathogenic_SNP_info, 2) .* ones(size(pathogenic_SNP_targets)));
pathogenic_error = (sum(sum(abs(out3 - pathogenic_SNP_targets) ./ (size(pathogenic_SNP_features, 2) .* ones(size(pathogenic_SNP_targets))))));

false_negative_group = find(out3 <= pathogenic_ratio);
true_positive_group  = find(out3 > pathogenic_ratio);

number_of_pathogenic_SNP = size(pathogenic_SNP_features, 2);

number_of_false_negative = size(false_negative_group, 2);
number_of_true_positive  = size(true_positive_group, 2);

%------------------------------- reserved ---------------------------------
%{
false_negative_pathogenic_SNP_info     = pathogenic_SNP_info(false_negative_group);
false_negative_pathogenic_SNP_features = pathogenic_SNP_features(:, false_negative_group);
false_negative_pathogenic_SNP_targets  = pathogenic_SNP_targets(:, false_negative_group);

true_positive_pathogenic_SNP_info     = pathogenic_SNP_info(true_positive_group);
true_positive_pathogenic_SNP_features = pathogenic_SNP_features(:, true_positive_group);
true_positive_pathogenic_SNP_targets  = pathogenic_SNP_targets(:, true_positive_group);

[out1, out2] = forward_probabilistic_pass_2layer(best_weight_first_layer, best_weight_second_layer, false_negative_pathogenic_SNP_features);
error = evaluate_probabilistic_model(out2, false_negative_pathogenic_SNP_targets)

[out1, out2] = forward_probabilistic_pass_2layer(best_weight_first_layer, best_weight_second_layer, true_positive_pathogenic_SNP_features);
error = evaluate_probabilistic_model(out2, true_positive_pathogenic_SNP_targets)
%}
%------------------------------- reserved ---------------------------------
accuracy          = (number_of_true_negative + number_of_true_positive) / (number_of_true_negative + number_of_true_positive + number_of_false_negative + number_of_false_positive);
sensitivity       = number_of_true_positive / (number_of_true_positive + number_of_false_negative);
specificity       = number_of_true_negative / (number_of_true_negative + number_of_false_positive);
balanced_accuracy = (sensitivity + specificity) / 2;

balanced_error    = (neutral_error + pathogenic_error) / 2;

disp(sprintf('neutral SNPs        : %d', number_of_neutral_SNP));
disp(sprintf('pathogenic SNPs     : %d', number_of_pathogenic_SNP));
disp(sprintf('pathogenic ratio    : %6.4f', pathogenic_ratio));
disp('***************************');
disp(sprintf('false positive SNPs : %d', number_of_false_positive));
disp(sprintf('true negative SNPs  : %d', number_of_true_negative));
disp(sprintf('false negative SNPs : %d', number_of_false_negative));
disp(sprintf('true positive SNPs  : %d', number_of_true_positive));
disp('***************************');
disp(sprintf('accuracy            : %6.4f', accuracy));
disp(sprintf('sensitivity         : %6.4f', sensitivity));
disp(sprintf('specificity         : %6.4f', specificity));
disp(sprintf('balanced accuracy   : %6.4f', balanced_accuracy));
disp('***************************');
disp(sprintf('neutral probabilistic error    : %6.4f', neutral_error));
disp(sprintf('pathogenic probabilistic error : %6.4f', pathogenic_error));
disp(sprintf('balanced probabilistic error   : %6.4f', balanced_error));

neutral_SNP_index    = find(testing_targets == 0);
neutral_SNP_features = testing_features(:, neutral_SNP_index);
[number_of_features number_of_neutral_samples] = size(neutral_SNP_features);

[neutral_out1, neutral_out2, neutral_out3] = forward_pass(best_weight_first_layer, best_weight_second_layer, best_weight_third_layer, neutral_SNP_features);

pathogenic_SNP_index    = find(testing_targets == 1);
pathogenic_SNP_features = testing_features(:, pathogenic_SNP_index);
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


%{
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
title('Polyphen2');

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
title('Mutation Taster');

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
title('GERP++');

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
title('Combiner');
%}

%{
subplot(1,1,1);
neutral_distribution = hist(neutral_out3, [0:0.02:1]);
plot([0:0.02:1], neutral_distribution, '-.b');
hold on
pathogenic_distribution = hist(pathogenic_out3, [0:0.02:1]);
plot([0:0.02:1], pathogenic_distribution, '-.r');
legend('neutral variants score', 'pathogenic variants score');
hold off
ylabel('samples');
xlabel('score');
title('Combiner');
%}

set(findall(f,'type','text'),'fontSize',14,'fontName','Times New Roman');
%set(findall(f,'type','text'),'fontWeight','bold');
%close(f);
