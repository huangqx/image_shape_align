function [Corres] = i2s_nn_query_pruned(...
    targetPC_2D,...
    sourcePC_2D,...
    minSigmaCorres,...
    lambda_saliency)
% the third dimension of the sourcePC_2D and targetPC_2D describe the
% saliency information of the point clouds (derived from the edgeDetect)
% argument 'lambda_saliency' is used to scale this coordinate
% argument 'sigmaCorres' is used to weight the nearest-neighbor
% correspondences

sourcePC_2D(3,:) = sourcePC_2D(3,:)*lambda_saliency;
targetPC_2D(3,:) = targetPC_2D(3,:)*lambda_saliency;

[footIds_s_in_t, d_s_in_t] = knnsearch(targetPC_2D', sourcePC_2D');
[footIds_t_in_s, d_t_in_s] = knnsearch(sourcePC_2D', targetPC_2D');

sigma = max(minSigmaCorres, median(d_s_in_t));
fprintf('sigma = %f.\n', sigma);
w_s_in_t = sigma^2./(d_s_in_t'.^2 + sigma^2);
w_t_in_s = sigma^2./(d_t_in_s'.^2 + sigma^2);

ids_s_in_t = find(d_s_in_t < 4*sigma);
ids_t_in_s = find(d_t_in_s < 4*sigma);
Corres = [ids_s_in_t', footIds_t_in_s(ids_t_in_s)';
    footIds_s_in_t(ids_s_in_t)', ids_t_in_s'];

Corres = [Corres; max(1e-2,[w_s_in_t(ids_s_in_t), w_t_in_s(ids_t_in_s)])];
saliencyWeights = sqrt(sourcePC_2D(3, Corres(1,:)).*...
    targetPC_2D(3, Corres(2,:)))/lambda_saliency;
Corres(3,:) = Corres(3,:).*saliencyWeights;
