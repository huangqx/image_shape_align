function [Camera_opt, partTrans_opt] = mm_part_affine_align(Shape_init,...
            partTrans_init,...
            partCorres,...
            Camera_init,...
            Para)
% This function optimizes an axis-aligned affine transformation for each
% part so that the deformed shape aligns with the image-correspondences
% Input arguments:
%   Shape_init: Input image
%   partTrans_init: the initial part affine transformations
%   partCorres: The correspondences between points on part and their
%               corresponding points in 2D
%   Camera_init: the initial camera configuration
%   Para: please refer to the software document
% Output arguments:
%   Camera_opt: the optimized camera configuration
%   partTrans_opt: the optimized part affine transformations

axis_z = Camera_init.origin - Camera_init.lookAt;
viewDis = norm(axis_z); % ViewDistance
axis_z = axis_z/viewDis;
axis_y = Camera_init.upVec;
axis_y = axis_y - axis_z*(axis_z'*axis_y);
axis_y = axis_y/norm(axis_y);
axis_x = cross(axis_y, axis_z);

% Write out the camera variables to be optimized
Var.s = Camera_init.scale;
Var.R = [axis_x, axis_y, axis_z]';
Var.o = Camera_init.lookAt;
Var.f = viewDis;

% Initialize the transformation parameters
Var.Rot_pairs = partTrans_init.Rot_pairs;
Var.Trans_pairs = partTrans_init.Trans_pairs;
Var.theta_sym = partTrans_init.theta_sym;
Var.Trans_yz = partTrans_init.Trans_yz;

numPairs = size(Var.Trans_pairs, 2);
numSyms = size(Var.theta_sym, 2);
numVars = 7 + 3*numPairs + 2*numSyms;

g = zeros(numVars, 1);
H = zeros(numVars, numVars);

for iter = 1:Para.numIterations_non_rigid
    % Generate the deformed source point cloud
    for i = 1:numSyms
        partId = Shape_init.nonsym_ids(i);
        pc = partCorres{partId};
        theta = Var.theta_sym(i);
        R = [1, 0, 0;
            0, cos(theta), -sin(theta);
            0, sin(theta), cos(theta)];
        t = [0, Var.Trans_yz(1, i), Var.Trans_yz(2,i)]';
        vertexPoss = Shape_init.vertexPoss(:, Shape_init.meshes{partId}.vertexIds);
        [g_p, H_p] = part_disp(pc, Var, R, t, vertexPoss, Para.lambda_1st_part);
        ids = [1:7, 9:10];
        g_p = g_p(ids);
        H_p = H_p(ids, ids);
        ids = [1:7, (2*i+6):(2*i+7)];
        H(ids,ids) = H(ids,ids) + H_p;
        g(ids) = g(ids) + g_p;
    end
    for i = 1:numPairs
        part1Id = Shape_init.sym_pairs(1, i);
        pc = partCorres{part1Id};
        R = Var.Rot_pairs(:,:,i);
        t = Var.Trans_pairs(:,i);
        vertexPoss = Shape_init.vertexPoss(:,...
            Shape_init.meshes{part1Id}.vertexIds);
        [g_p, H_p] = part_disp(pc, Var, R, t, vertexPoss, Para.lambda_1st_part);
        g_p = g_p(1:10);
        H_p = H_p(1:10,1:10);
        ids = [1:7, (2*numSyms+3*i+5):(2*numSyms+3*i+7)];
        H(ids,ids) = H(ids,ids) + H_p;
        g(ids) = g(ids) + g_p;
 
        part2Id = Shape_init.sym_pairs(2, i);
        pc = partCorres{part2Id};
        
        vertexPoss = Shape_init.vertexPoss(:,...
            Shape_init.meshes{part2Id}.vertexIds);
        
        t(1) = -t(1);
        R = diag([-1,1,1])*R*diag([-1,1,1]);
        [g_p, H_p] = part_disp(pc, Var, R, t, vertexPoss, Para.lambda_1st_part);
        D = diag([ones(1,7), [-1,1,1,1,-1,-1]]);
        H_p = D*H_p*D;
        g_p = D*g_p;
        g_p = g_p(1:10);
        H_p = H_p(1:10,1:10);
        H(ids,ids) = H(ids,ids) + H_p;
        g(ids) = g(ids) + g_p;
    end
    
    x = H\g;
    
    % Perform line search to find the optimal solution
    % Calculate the bejective function at the previous iteration
    energy = obj_val(partCorres, Var, Shape_init, Para.lambda_1st_part);
    
    alpha = 1;
    succces = 0;
    for searchIter = 1:10
        Var1.s = Var.s;
        Var1.f = Var.f + x(1)*alpha;
        Var1.o = Var.o + x(5:7)*alpha;
        C = [0, -x(4), x(3);
            x(4), 0, -x(2);
            -x(3), x(2), 0]*alpha;
        Var1.R = normalize_rot(expm(C)*Var.R);
        
        for i = 1:numSyms
            dx = x((2*i+6):(2*i+7))*alpha;
            Var1.Trans_yz(:,i) = Var.Trans_yz(:,i) + dx(1:2);
            Var1.theta_sym(i) = 0 + Var.theta_sym(i);
        end
        for i = 1:numPairs
            dx = x((2*numSyms + 3*i+5):(2*numSyms + 3*i+7))*alpha;
%            dR = expm([0, -dx(6), dx(5);
%                      dx(6), 0, -dx(4);
%                      -dx(5), dx(4), 0]);
%            Var1.Rot_pairs(:,:,i) = normalize_rot(dR*Var.Rot_pairs(:,:,i));
            Var1.Rot_pairs(:,:,i) = Var.Rot_pairs(:,:,i);
            Var1.Trans_pairs(:,i) = Var.Trans_pairs(:,i) + dx(1:3); 
        end
    
        % Calculate the distance after displacing the shapes
        e = obj_val(partCorres, Var1, Shape_init, Para.lambda_1st_part);
        if e < energy
            Var = Var1;
            succces = 1;
            break;
        end
        alpha = alpha/2;
    end
    if succces == 0
        break;
    end
end

% Write back the optimized camera configuration
Camera_opt = Camera_init;
Camera_opt.scale = Var.s;
Camera_opt.lookAt = Var.o;
Camera_opt.origin = Var.o + Var.R(3,:)'*Var.f;
Camera_opt.upVec = Var.R(2,:)';

partTrans_opt.Rot_pairs = Var.Rot_pairs;
partTrans_opt.Trans_pairs = Var.Trans_pairs;
partTrans_opt.theta_sym = Var.theta_sym;
partTrans_opt.Trans_yz = Var.Trans_yz;

function [e] = obj_val(partCorres, Var, Shape_init, lambda_1st)
e = 0;
numSyms = length(Shape_init.nonsym_ids);
numPairs = size(Shape_init.sym_pairs, 2);
for i = 1:numSyms
    partId = Shape_init.nonsym_ids(i);
    pc = partCorres{partId};
    theta = Var.theta_sym(i);
    R = [1, 0, 0;
        0, cos(theta), -sin(theta);
        0, sin(theta), cos(theta)];
    t = [0, Var.Trans_yz(1, i), Var.Trans_yz(2,i)]';
    vertexPoss = Shape_init.vertexPoss(:, Shape_init.meshes{partId}.vertexIds);
    e_part = obj_val_part(pc, Var, R, t, vertexPoss, lambda_1st);
    e = e + e_part;
end
for i = 1:numPairs
    part1Id = Shape_init.sym_pairs(1, i);
    pc = partCorres{part1Id};
    R = Var.Rot_pairs(:,:,i);
    t = Var.Trans_pairs(:,i);
    vertexPoss = Shape_init.vertexPoss(:,...
        Shape_init.meshes{part1Id}.vertexIds);
    e_part = obj_val_part(pc, Var, R, t, vertexPoss, lambda_1st);
    e = e + e_part;
 
    part2Id = Shape_init.sym_pairs(2, i);
    pc = partCorres{part2Id};
        
    vertexPoss = Shape_init.vertexPoss(:,...
        Shape_init.meshes{part2Id}.vertexIds);
        
    t(1) = -t(1);
    R = diag([-1,1,1])*R*diag([-1,1,1]);
    e_part = obj_val_part(pc, Var, R, t, vertexPoss, lambda_1st);
    e = e + e_part;
end
    


function [e] = obj_val_part(pc, Var, R_cur, t_cur, vertexPoss, lambda_1st)
%
lambda = lambda_1st*size(pc.points3d, 2)/size(vertexPoss, 2);

numV = size(vertexPoss, 2);
vertexPoss_cur = R_cur*vertexPoss + t_cur*ones(1, numV);

% The regularization term
d = vertexPoss_cur - vertexPoss;
e = lambda*sum(sum(d.*d));

% The data-term
numP = size(pc.points3d, 2);
sourcePC_3D = R_cur*pc.points3d + t_cur*ones(1, numP);

P_bar = Var.R*(sourcePC_3D - Var.o*ones(1,numP));
P_hat(1,:) = Var.f*P_bar(1,:)./(Var.f - P_bar(3,:));
P_hat(2,:) = Var.f*P_bar(2,:)./(Var.f - P_bar(3,:));
P_hat = P_hat/Var.s;

w = kron(pc.weights, ones(1,2));
W = sparse(1:length(w), 1:length(w), w);

r = zeros(2*numP, 1);

r(1:2:(2*numP)) = P_hat(1,:)'- pc.points2d(1,:)';
r(2:2:(2*numP)) = P_hat(2,:)'- pc.points2d(2,:)';
e = e + r'*W*r;

function [g, H] = part_disp(pc, Var, R_cur, t_cur, vertexPoss, lambda_1st)
%
H = zeros(13, 13);
g = zeros(13, 1);

lambda = lambda_1st*size(pc.points3d, 2)/size(vertexPoss, 2);

numV = size(vertexPoss, 2);
vertexPoss_cur = R_cur*vertexPoss + t_cur*ones(1, numV);

% The regularization term
Jr = kron(ones(numV,1), [eye(3), zeros(3,3)]);
gr = reshape(vertexPoss - vertexPoss_cur, [3*numV,1]);

for i = 1:size(vertexPoss, 2)
    Jr(1:3:(3*numV), 5) = vertexPoss_cur(3,:)';
    Jr(1:3:(3*numV), 6) = -vertexPoss_cur(2,:)';
    Jr(2:3:(3*numV), 4) = -vertexPoss_cur(3,:)';
    Jr(2:3:(3*numV), 6) = vertexPoss_cur(1,:)';
    Jr(3:3:(3*numV), 4) = vertexPoss_cur(2,:)';
    Jr(3:3:(3*numV), 5) = -vertexPoss_cur(1,:)';
end
H(8:13,8:13) = lambda*(Jr'*Jr);
g(8:13) = lambda*Jr'*gr;

% The data-term
numP = size(pc.points3d, 2);
sourcePC_3D = R_cur*pc.points3d + t_cur*ones(1, numP);

P_bar = Var.R*(sourcePC_3D - Var.o*ones(1,numP));
P_hat(1,:) = Var.f*P_bar(1,:)./(Var.f - P_bar(3,:));
P_hat(2,:) = Var.f*P_bar(2,:)./(Var.f - P_bar(3,:));
P_hat = P_hat/Var.s;

w = kron(pc.weights, ones(1,2));
W = sparse(1:length(w), 1:length(w), w);

J2 = zeros(3*numP, 6);
J = zeros(2*numP, 14);
J3 = kron(ones(numP, 1), [eye(3), zeros(3,3)]);
r = zeros(2*numP, 1);

for i = 1:size(vertexPoss, 2)
    J3(1:3:(3*numP), 5) = P_bar(3,:)';
    J3(1:3:(3*numP), 6) = -P_bar(2,:)';
    J3(2:3:(3*numP), 4) = -P_bar(3,:)';
    J3(2:3:(3*numP), 6) = P_bar(1,:)';
    J3(3:3:(3*numP), 4) = P_bar(2,:)';
    J3(3:3:(3*numP), 5) = -P_bar(1,:)';
end

J3 = kron(sparse(1:numP, 1:numP, ones(1,numP)), Var.R)*J3;

% 
J2(:,4:6) = -kron(ones(numP,1), Var.R);
J2(1:3:(3*numP), 2) = P_bar(3,:)';     J2(1:3:(3*numP), 3) = -P_bar(2,:)';
J2(2:3:(3*numP), 1) = -P_bar(3,:)';    J2(2:3:(3*numP), 3) = P_bar(1,:)';
J2(3:3:(3*numP), 1) = P_bar(2,:)';     J2(3:3:(3*numP), 2) = -P_bar(1,:)';
r(1:2:(2*numP)) = P_hat(1,:)'- pc.points2d(1,:)';
r(2:2:(2*numP)) = P_hat(2,:)'- pc.points2d(2,:)';
%
J(1:2:(2*numP),1) = ((Var.f*P_bar(1,:))./(Var.f-P_bar(3,:)))'/(Var.s^2);
J(2:2:(2*numP),1) = ((Var.f*P_bar(2,:))./(Var.f-P_bar(3,:)))'/(Var.s^2);
%
J(1:2:(2*numP),2) = ((P_bar(1,:).*P_bar(3,:))./((Var.f-P_bar(3,:)).^2))'/Var.s;
J(2:2:(2*numP),2) = ((P_bar(2,:).*P_bar(3,:))./((Var.f-P_bar(3,:)).^2))'/Var.s;
%
c1 = -kron((Var.f./(Var.f - P_bar(3,:)'))/Var.s, ones(1,6));
c2 = -kron(((Var.f*P_bar(1,:)')./((Var.f - P_bar(3,:)').^2))/Var.s,ones(1,6));
J(1:2:(2*numP),3:8) = J2(1:3:(3*numP),:).*c1 + J2(3:3:(3*numP),:).*c2;
c1 = -kron((Var.f./(Var.f - P_bar(3,:)'))/Var.s, ones(1,6));
c2 = -kron(((Var.f*P_bar(1,:)')./((Var.f - P_bar(3,:)').^2))/Var.s,ones(1,6));
J(1:2:(2*numP),9:14) = c1.*J3(1:3:(3*numP),:) +c2.*J3(3:3:(3*numP),:);
        
c1 = -kron((Var.f./(Var.f - P_bar(3,:)'))/Var.s, ones(1,6));
c2 = -kron(((Var.f*P_bar(2,:)')./((Var.f - P_bar(3,:)').^2))/Var.s,ones(1,6));
J(2:2:(2*numP),3:8) = J2(2:3:(3*numP),:).*c1 + J2(3:3:(3*numP),:).*c2;
c1 = -kron((Var.f./(Var.f - P_bar(3,:)'))/Var.s, ones(1,6));
c2 = -kron(((Var.f*P_bar(2,:)')./((Var.f - P_bar(3,:)').^2))/Var.s,ones(1,6));
J(2:2:(2*numP),9:14) = c1.*J3(2:3:(3*numP),:) +c2.*J3(3:3:(3*numP),:);
    
% Generate the Gaussian-Newton approximant of the quadratic energy
J = J(:, 2:14);
H = H + J'*W*J;
g = g + J'*W*r;


function [R] = normalize_rot(R_in)
R_in(:,1) = R_in(:,1)/norm(R_in(:,1));
R_in(:,2) = R_in(:,2) - R_in(:,1)*(R_in(:,2)'*R_in(:,1));
R_in(:,2) = R_in(:,2)/norm(R_in(:,2));
R_in(:,3) = R_in(:,3) - R_in(:,1)*(R_in(:,3)'*R_in(:,1)) - R_in(:,2)*(R_in(:,3)'*R_in(:,2));
R_in(:,3) = R_in(:,3)/norm(R_in(:,3));
R = R_in;
