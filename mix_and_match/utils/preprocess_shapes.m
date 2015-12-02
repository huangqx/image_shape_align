function [Shapes_out] = preprocess_shapes(Shapes)
% Detect part structures of shapes
% Align them in a common coordinate system
%

for i = 1:length(Shapes)
    % normalize the shape
    Shapes{i} = normalize_shape(Shapes{i});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aligning stage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default aligning parameters
Para_align.gridRes = 0.1500;
Para_align.numSamples = 8192;
Para_align.lambda_first = 0.2500;
Para_align.lambda_smooth = 2.5000;
Para_align.numIterations_outer = 8;
Para_align.numIterations_alternate = 16;
Para_align.numIterations_pairwise = 25;

[SAMPLE, PAIRMATCH] = all_pairwise_align(Shapes, Para_align);
FFDs_opt = joint_align(Shapes, SAMPLE, PAIRMATCH, Para_align);

for i = 1 : length(Shapes)
    C_source_ori = sp_ffd_basis_coeff(FFDs_opt{i}, Shapes{i}.vertexPoss);
    % FFD coefficient
    Shapes{i}.vertexPoss = FFDs_opt{i}.ctrlPos_cur*C_source_ori';
end
Shapes_out = Shapes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Shape_out] = normalize_shape(Shape)
%
boxMin = min(Shape.vertexPoss')';
boxMax = max(Shape.vertexPoss')';
center = (boxMin + boxMax)/2;
scale = 1/norm(boxMax - boxMin);
numVertices = size(Shape.vertexPoss, 2);
Shape.vertexPoss = Shape.vertexPoss - center*ones(1, numVertices);
Shape.vertexPoss = Shape.vertexPoss*scale;
Shape_out = Shape;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SAMPLE, PAIRMATCH] = all_pairwise_align(Shapes,... % The input shapes
    Para_align) % The parameters
% This function performs all pair-wise matching along pairs of shapes
% specified by a shape graph

% Perform sampling on the input shapes
numShapes = length(Shapes);
for id = 1 : numShapes
    SAMPLE{id} = sp_mesh_sampling(Shapes{id}, Para_align.numSamples);
end

for id = 1 : numShapes
    FFDs{id} = sp_ffd_init_sym(Shapes{id}, Para_align.gridRes);
end

edges = zeros(2, numShapes*(numShapes-1)/2);
count = 0;
for sId = 1 : numShapes
    for tId = (sId+1): numShapes
        count = count + 1;
        edges(1, count) = sId;
        edges(2, count) = tId;
    end
end

% Perform pair-wise alignment
for pairId = 1 : size(edges, 2)
    PAIRMATCH{pairId} = pairwise_matching(...
        FFDs{edges(1, pairId)},...
        SAMPLE{edges(1, pairId)},...
        SAMPLE{edges(2, pairId)},...
        edges(1, pairId),...
        edges(2, pairId),...
        Para_align);
end

function [match] = pairwise_matching(...
    sourceFFD,...
    sourceSample,...
    targetSample,...
    sId, tId,...
    Para_align)
%
[ffd_opt, medianDis] = pairwise_ffd_align(...
    sourceFFD,...
    sourceSample,...
    targetSample,...
    Para_align);
    
% After alignment, compute the closest point pairs
C_source_ori = sp_ffd_basis_coeff(ffd_opt, sourceSample); % FFD coefficient
sourcePoints = ffd_opt.ctrlPos_cur*C_source_ori'; % Current point positions
[IDX_s_t, dis_s_t] = knnsearch(targetSample', sourcePoints'); % Compute nearest neighbors
sigma = median(dis_s_t); % Use robust thresholding to find correspondences
ids = find(dis_s_t < sigma*2);
corres = [ids'; IDX_s_t(ids)'];
    
% Store the correspondences
match.sId = sId;
match.tId = tId;
match.corres = corres;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pair-wise ffd alignment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ffd_opt, medianDis] = pairwise_ffd_align(...
    sourceFFD, sourceSample, targetSample, Para)
% Optimize a free-form deformation that aligns the source shape to the
% target shape

% The FFD coefficient of the source points
C_source_ori = sp_ffd_basis_coeff(sourceFFD, sourceSample);

% Applies non-rigid 
for iter = 1:Para.numIterations_pairwise
    sourcePoints = sourceFFD.ctrlPos_cur*C_source_ori';
    % Compute correspondences
    [Corres, medianDis] = sp_closest_point_corres(sourcePoints,...
        targetSample);
    
    % Deform the shape accordingly
    nC = size(Corres, 2);
    W = sparse(1:nC, 1:nC, Corres(3,:));
    
    Ds = C_source_ori(Corres(1,:),:)';
    Pt = targetSample(:, Corres(2,:));
    
    dimX = size(Ds, 1);
    A = Ds*W*Ds' + Para.lambda_first*eye(dimX)...
        + Para.lambda_smooth*sourceFFD.H_smooth;
    b = Ds*W*Pt' + Para.lambda_first*sourceFFD.ctrlPos_ori';
    sourceFFD.ctrlPos_cur = (A\b)';
end
%
ffd_opt = sourceFFD;

%C_source_ori = sp_ffd_basis_coeff(sourceFFD, sourceSample);
%sourceSample_def = sourceFFD.ctrlPos_cur*C_source_ori';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute closest point pairs between three point clouds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Corres, medianDis] = sp_closest_point_corres(sourcePoints,...
    targetPoints)
%
% Compute cloest point pairs in both directions
[footIds_s_in_t, d_s_in_t] = knnsearch(targetPoints', sourcePoints');
[footIds_t_in_s, d_t_in_s] = knnsearch(sourcePoints', targetPoints');

dis = [d_s_in_t', d_t_in_s'];
Corres_s = [1:length(footIds_s_in_t), footIds_t_in_s'];
Corres_t = [footIds_s_in_t', 1:length(footIds_t_in_s)];

Corres = [Corres_s; Corres_t];
d = sort(dis);
sigma = d(floor(length(d)*0.85));
sigma = max(5e-2, sigma);
weights = sigma*sigma./(sigma*sigma + dis.*dis);

Corres = [Corres; weights];
medianDis = sqrt(mean(dis.*dis));

function [FFDs] = joint_align(Shapes, SAMPLE, PAIRMATCH, Para_align)
% Given the correspondences computed between all pairs of shapes
% Align all the shapes in the world coordinate system
% Input arguments:
%       Shapes: the input shapes
%       SAMPLE: the samples placed on each shape
%       PAIRMATCH: the pre-computed correspondences between pairs of shapes
%       Para_align: the parameters
% Output argument:
%        FFDs{shapeId}: the optimized free-form deformation of each shape

numShapes = length(SAMPLE);
numPairs = length(PAIRMATCH);

% Allocate space to store the pair-wise correspondence weights
for pairId = 1 : numPairs
    nc = size(PAIRMATCH{pairId}.corres, 2);
    PAIRMATCH{pairId}.weights = ones(1, nc);
end

% Count the number of total correspondences
numC_all = 0;
for pairId = 1 : numPairs
    numC = size(PAIRMATCH{pairId}.corres, 2);
    numC_all = numC_all + numC;
end

for outIter = 1:Para_align.numIterations_outer
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Precompute the number of correspondences acted on each sample
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Allocate space for cumulative weights
    for shapeId = 1 : numShapes
        weights{shapeId} = zeros(1, size(SAMPLE{shapeId}, 2));
    end

    % Pre-compute the weights
    for pairId = 1 : numPairs
        nc = size(PAIRMATCH{pairId}.corres, 2);
        PAIRMATCH{pairId}.midPoints = zeros(3, nc);
        sId = PAIRMATCH{pairId}.sId;
        tId = PAIRMATCH{pairId}.tId;
        ids1 = PAIRMATCH{pairId}.corres(1,:);
        ids2 = PAIRMATCH{pairId}.corres(2,:);
        weights{sId}(ids1) = weights{sId}(ids1) + PAIRMATCH{pairId}.weights;
        weights{tId}(ids2) = weights{tId}(ids2) + PAIRMATCH{pairId}.weights;
    end

    % Pre-compute the deformation structure used for each shape
    SAMPLE_def = SAMPLE;
    for shapeId = 1 : numShapes
        FFDs{shapeId} = sp_ffd_init_sym(Shapes{shapeId}, Para_align.gridRes);
        Term{shapeId}.b = sp_ffd_basis_coeff(FFDs{shapeId}, SAMPLE{shapeId});
        w = weights{shapeId};
        W = sparse(1:length(w), 1:length(w), w);
        dimX = size(Term{shapeId}.b, 2);
        Term{shapeId}.A = Term{shapeId}.b'*W*Term{shapeId}.b...
            + Para_align.lambda_first*eye(dimX)...
            + Para_align.lambda_smooth*FFDs{shapeId}.H_smooth;
        Term{shapeId}.b = Term{shapeId}.b';
    end
    
    % Peform alternating optimization to optimize the deformation on each
    % shape
    for iter = 1:Para_align.numIterations_alternate
        % Compute the deformed positions of sample points
        for shapeId = 1 : numShapes
            SAMPLE_def{shapeId} = FFDs{shapeId}.ctrlPos_cur*Term{shapeId}.b;
        end
    
        % Compute the intermediate points
        sqrDis = 0;
        for pairId = 1 : numPairs
            sId = PAIRMATCH{pairId}.sId;
            tId = PAIRMATCH{pairId}.tId;
            poss_s = SAMPLE_def{sId}(:, PAIRMATCH{pairId}.corres(1,:));
            poss_t = SAMPLE_def{tId}(:, PAIRMATCH{pairId}.corres(2,:));
            PAIRMATCH{pairId}.midPoints = (poss_s + poss_t)/2;
            d = poss_s - poss_t;
            sqrDis = sqrDis + sum(sum(d.*d).*PAIRMATCH{pairId}.weights);
        end
 %       fprintf('sqrDis = %f.\n', sqrDis);
    
        % Perform the alignment
        for shapeId = 1 : numShapes
            TP{shapeId} = zeros(3, size(SAMPLE{shapeId}, 2));
        end
        for pairId = 1 : numPairs
            sId = PAIRMATCH{pairId}.sId;
            tId = PAIRMATCH{pairId}.tId;
            ids1 = PAIRMATCH{pairId}.corres(1,:);
            ids2 = PAIRMATCH{pairId}.corres(2,:);
            buf = PAIRMATCH{pairId}.midPoints.*(ones(3,1)*PAIRMATCH{pairId}.weights);
            TP{sId}(:, ids1) = TP{sId}(:, ids1) + buf;
            TP{tId}(:, ids2) = TP{tId}(:, ids2) + buf;
        end
        for shapeId = 1 : numShapes
            b = Term{shapeId}.b*TP{shapeId}' +...
                Para_align.lambda_first*FFDs{shapeId}.ctrlPos_ori';
            FFDs{shapeId}.ctrlPos_cur = (Term{shapeId}.A\b)';
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reweight the correspondences
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for shapeId = 1 : numShapes
        SAMPLE_def{shapeId} = FFDs{shapeId}.ctrlPos_cur*Term{shapeId}.b;
    end
    %
    sqrDisVec = zeros(1, numC_all);
    off = 0;
    for pairId = 1 : numPairs
        sId = PAIRMATCH{pairId}.sId;
        tId = PAIRMATCH{pairId}.tId;
        poss_s = SAMPLE_def{sId}(:, PAIRMATCH{pairId}.corres(1,:));
        poss_t = SAMPLE_def{tId}(:, PAIRMATCH{pairId}.corres(2,:));
        d = poss_s - poss_t;
        d = sum(d.*d);
        nc = length(d);
        sqrDisVec((off+1):(off+nc)) = d;
        off = off + nc;
    end
    sigma = median(sqrDisVec);
    for pairId = 1 : numPairs
        sId = PAIRMATCH{pairId}.sId;
        tId = PAIRMATCH{pairId}.tId;
        poss_s = SAMPLE_def{sId}(:, PAIRMATCH{pairId}.corres(1,:));
        poss_t = SAMPLE_def{tId}(:, PAIRMATCH{pairId}.corres(2,:));
        d = poss_s - poss_t;
        d = sum(d.*d);
        PAIRMATCH{pairId}.weights = sqrt(sigma)./sqrt(sigma + d);
    end
end

