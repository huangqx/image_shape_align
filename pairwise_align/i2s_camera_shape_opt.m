function [Camera_opt, ctrlPos_opt] = i2s_camera_shape_opt(...
        targetPC_2D,... % The 2D edge points extracted from the input image
        sourcePC_3D_ori,... % The corresponding 3D points extracted from the input shape
        corresWeights,... % Weights associated with the correspondences
        Shape_FFD,...  % The Symmetric FFD deformation structure assocaited with the input shape
        Camera_init,... % The initial camera configuration
        Para)
% Para.numIterations_non_rigid: the number of gauss-newton iterations
%                               (default value is 10)
% Para.lambda_1st: penalize the absoluate deformation
% Para.lambda_smooth: make sure that the deformation is smooth
    
% Compute the camera coordinate system
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

% The deformation field to be optimized
Var.ctrlPos_opt = Shape_FFD.ctrlPos_cur;


% Pre-computed ffd basis coefficient for the 3D-points
numP = size(targetPC_2D, 2);
ffd_coeff_sourcePC = i2s_ffd_basis_coeff(Shape_FFD,...
    double(sourcePC_3D_ori));
numCtrl = size(Var.ctrlPos_opt, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The constraint set, which enforces that the deformation is symmetric
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A, b] = combine_constraints(Shape_FFD);
A = [sparse(size(A,1), 7), A];


% Allocate variables to hold the residual and the jacobi
P_hat = zeros(2, numP);
J2 = zeros(3*numP, 6); % A buffer for the computing the jacobi matrix
J = zeros(2*numP, 8+3*numCtrl);  % store the jacobi matrix
r = zeros(2*numP, 1); % store the residual of correspondences

% Convert the weight vector into a diagonal matrix
W = sparse(1:(2*numP),1:(2*numP), kron(corresWeights, ones(1,2)));


for iter = 1:Para.numIterations_non_rigid
    % Generate the deformed source point cloud
    sourcePC_3D = Var.ctrlPos_opt*ffd_coeff_sourcePC';
   
    % Apply the camera parameter to obtain the 2D location
    P_bar = Var.R*(sourcePC_3D - Var.o*ones(1,numP));
    P_hat(1,:) = Var.f*P_bar(1,:)./(Var.f - P_bar(3,:));
    P_hat(2,:) = Var.f*P_bar(2,:)./(Var.f - P_bar(3,:));
    P_hat = P_hat/Var.s;
    
    % 
    J2(:,4:6) = -kron(ones(numP,1), Var.R);
    J2(1:3:(3*numP), 2) = P_bar(3,:)';     J2(1:3:(3*numP), 3) = -P_bar(2,:)';
    J2(2:3:(3*numP), 1) = -P_bar(3,:)';    J2(2:3:(3*numP), 3) = P_bar(1,:)';
    J2(3:3:(3*numP), 1) = P_bar(2,:)';     J2(3:3:(3*numP), 2) = -P_bar(1,:)';
    r(1:2:(2*numP)) = P_hat(1,:)'- targetPC_2D(1,:)';
    r(2:2:(2*numP)) = P_hat(2,:)'- targetPC_2D(2,:)';
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
    c1 = -kron((Var.f./(Var.f - P_bar(3,:)'))/Var.s, ones(1,numCtrl));
    c2 = -kron(((Var.f*P_bar(1,:)')./((Var.f - P_bar(3,:)').^2))/Var.s,ones(1,numCtrl));
    J(1:2:(2*numP),9:(8+numCtrl)) = (c1*Var.R(1,1)+c2*Var.R(3,1)).*ffd_coeff_sourcePC;
    J(1:2:(2*numP),(9+numCtrl):(8+2*numCtrl)) = (c1*Var.R(1,2)+c2*Var.R(3,2)).*ffd_coeff_sourcePC;
    J(1:2:(2*numP),(9+2*numCtrl):(8+3*numCtrl)) = (c1*Var.R(1,3)+c2*Var.R(3,3)).*ffd_coeff_sourcePC;
    
    
    c1 = -kron((Var.f./(Var.f - P_bar(3,:)'))/Var.s, ones(1,6));
    c2 = -kron(((Var.f*P_bar(2,:)')./((Var.f - P_bar(3,:)').^2))/Var.s,ones(1,6));
    J(2:2:(2*numP),3:8) = J2(2:3:(3*numP),:).*c1 + J2(3:3:(3*numP),:).*c2;
    c1 = -kron((Var.f./(Var.f - P_bar(3,:)'))/Var.s, ones(1,numCtrl));
    c2 = -kron(((Var.f*P_bar(2,:)')./((Var.f - P_bar(3,:)').^2))/Var.s,ones(1,numCtrl));
    J(2:2:(2*numP),9:(8+numCtrl)) = (c1*Var.R(2,1)+c2*Var.R(3,1)).*ffd_coeff_sourcePC;
    J(2:2:(2*numP),(9+numCtrl):(8+2*numCtrl)) = (c1*Var.R(2,2)+c2*Var.R(3,2)).*ffd_coeff_sourcePC;
    J(2:2:(2*numP),(9+2*numCtrl):(8+3*numCtrl)) = (c1*Var.R(2,3)+c2*Var.R(3,3)).*ffd_coeff_sourcePC;
    
    % Generate the Gaussian-Newton approximant of the quadratic energy
    H = J'*W*J;
    g = J'*W*r;
    
%    fprintf(' object_val = %f.\n', r'*W*r);
    % Add the smoothness term to the formulation
    [H_regu, g_regu] = deform_regu_term(Shape_FFD.H_smooth,...
        Shape_FFD.ctrlPos_ori,...
        Var.ctrlPos_opt,...
        Para.lambda_1st,...
        Para.lambda_smooth);
    
    H = H + H_regu;
    g = g + g_regu;
    
    % The first variable is always fixed
    H = H(2:(8+3*numCtrl),2:(8+3*numCtrl));
    g = g(2:(8+3*numCtrl));
    
    % vector b should be shifted since we are computing an offset
    vec_P = reshape(Var.ctrlPos_opt', [3*numCtrl,1]);
    b_shift = b - A(:,8:(7+3*numCtrl))*vec_P;
    
    % Perform linearly constrained optimization to obtain the optimal
    % displacement
    x = quad_prog(H, g, A, b_shift);
    
    % Perform line search to find the optimal solution
    % Calculate the bejective function at the previous iteration
    energy = r'*W*r + obj_val_regu(Shape_FFD.H_smooth,...
        Shape_FFD.ctrlPos_ori,...
        Var.ctrlPos_opt,...
        Para.lambda_1st,...
        Para.lambda_smooth);
    
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
        dP = reshape(x(8:(3*numCtrl+7)), [numCtrl, 3])';
        Var1.ctrlPos_opt = Var.ctrlPos_opt + dP*alpha;
    
        % Calculate the distance after displacing the shapes
        e = obj_val(targetPC_2D, W, Var1, ffd_coeff_sourcePC);
        e = e + obj_val_regu(Shape_FFD.H_smooth,...
            Shape_FFD.ctrlPos_ori,...
            Var1.ctrlPos_opt,...
            Para.lambda_1st,...
            Para.lambda_smooth);
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
ctrlPos_opt = Var.ctrlPos_opt;

d_ori = norm(Camera_opt.lookAt-Camera_init.lookAt);
d_eye = norm(Camera_opt.origin-Camera_init.origin);
d_upVec = norm(Camera_opt.upVec-Camera_init.upVec);
fprintf('d_ori = %f, d_eye = %f, d_upVec = %f.\n', d_ori, d_eye, d_upVec); 

function [R] = normalize_rot(R_in)
R_in(:,1) = R_in(:,1)/norm(R_in(:,1));
R_in(:,2) = R_in(:,2) - R_in(:,1)*(R_in(:,2)'*R_in(:,1));
R_in(:,2) = R_in(:,2)/norm(R_in(:,2));
R_in(:,3) = R_in(:,3) - R_in(:,1)*(R_in(:,3)'*R_in(:,1)) - R_in(:,2)*(R_in(:,3)'*R_in(:,2));
R_in(:,3) = R_in(:,3)/norm(R_in(:,3));
R = R_in;

function [B, d] = combine_constraints(Shape_FFD)
%
[m1,n1] = size(Shape_FFD.Ax);
[m2,n2] = size(Shape_FFD.Ay);
[m3,n3] = size(Shape_FFD.Az);
m = m1 + m2 + m3;
n = n1 + n2 + n3;

B = sparse(m, n);
d = [Shape_FFD.bx; Shape_FFD.by; Shape_FFD.bz];
B(1:m1, 1:n1) = Shape_FFD.Ax;
B(m1+(1:m2), n1+(1:n2)) = Shape_FFD.Ay;
B(m1+m2+(1:m3), n1+n2+(1:n3)) = Shape_FFD.Az;

function [H_regu, g_regu] = deform_regu_term(...
    H_smooth, ctrlPos_ori, ctrlPos_cur,...
    lambda_1, lambda_2)

nP = size(ctrlPos_ori, 2);
H = H_smooth*lambda_2 + eye(nP)*lambda_1;
g = -lambda_2*H_smooth*ctrlPos_cur' +lambda_1*(ctrlPos_ori - ctrlPos_cur)';
H_regu = zeros(3*nP+8,3*nP+8);
g_regu = zeros(3*nP+8,1);
H_regu(9:(8+3*nP),9:(8+3*nP)) = kron(eye(3), full(H));
g_regu(9:(8+3*nP)) = reshape(g, [3*nP,1]);

function [e_regu] = obj_val_regu(...
    H_smooth, ctrlPos_ori, ctrlPos_cur,...
    lambda_1, lambda_2)

e1 = sum(sum((ctrlPos_ori-ctrlPos_cur).*(ctrlPos_ori-ctrlPos_cur)));
e2 = sum(diag(ctrlPos_cur*H_smooth*ctrlPos_cur'));
e_regu = lambda_1*e1 + lambda_2*e2;

function [x] = quad_prog(H,g, A, b)
% solve the following linearly constrained optimization problem
% min x'*H*x - 2*x'*g
% subject to A*x == b

nr = size(A, 1);
TP = [H, A'; A, zeros(nr, nr)];
t = [g;b];
x = TP\t;
x = x(1:size(H,1));

function [e] = obj_val(targetPC_2D, W, Var, ffd_coeff_sourcePC)

sourcePC_3D = Var.ctrlPos_opt*ffd_coeff_sourcePC';
numP = size(sourcePC_3D, 2);

% Apply the camera parameter to obtain the 2D location
P_bar = Var.R*(sourcePC_3D - Var.o*ones(1,numP));
P_hat(1,:) = Var.f*P_bar(1,:)./(Var.f - P_bar(3,:));
P_hat(2,:) = Var.f*P_bar(2,:)./(Var.f - P_bar(3,:));
P_hat = P_hat/Var.s;
    
%
r = zeros(2*numP, 1);
r(1:2:(2*numP)) = P_hat(1,:)'- targetPC_2D(1,:)';
r(2:2:(2*numP)) = P_hat(2,:)'- targetPC_2D(2,:)';
e = r'*W*r;
