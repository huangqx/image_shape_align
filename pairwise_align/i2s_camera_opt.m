function [Camera_opt] = i2s_camera_opt(...
    targetPC_2D,... % the feature edge points extracted in the natural object
    sourcePC_3D,... % the corresponding 3D points
    weights,... % the weights of the correspondences
    Camera_init,...
    Para) 

axis_z = Camera_init.origin - Camera_init.lookAt;
viewDis = norm(axis_z); % ViewDistance
axis_z = axis_z/viewDis;
axis_y = Camera_init.upVec;
axis_y = axis_y - axis_z*(axis_z'*axis_y);
axis_y = axis_y/norm(axis_y);
axis_x = cross(axis_y, axis_z);
%
Var.s = Camera_init.scale;
Var.R = [axis_x, axis_y, axis_z]';
Var.o = Camera_init.lookAt;
Var.f = viewDis;

numP = size(targetPC_2D, 2);

P_hat = zeros(2, numP);
J2 = zeros(3*numP, 6);
J = zeros(2*numP, 8);
g = zeros(2*numP, 1);

W = sparse(1:(2*numP),1:(2*numP), kron(weights, ones(1,2)));

for iter = 1:Para.numIterations_rigid
    P_bar = Var.R*(sourcePC_3D - Var.o*ones(1,numP));
    P_hat(1,:) = Var.f*P_bar(1,:)./(Var.f - P_bar(3,:));
    P_hat(2,:) = Var.f*P_bar(2,:)./(Var.f - P_bar(3,:));
    P_hat = P_hat/Var.s;
    J2(:,4:6) = -kron(ones(numP,1), Var.R);
    J2(1:3:(3*numP), 2) = P_bar(3,:)';     J2(1:3:(3*numP), 3) = -P_bar(2,:)';
    J2(2:3:(3*numP), 1) = -P_bar(3,:)';    J2(2:3:(3*numP), 3) = P_bar(1,:)';
    J2(3:3:(3*numP), 1) = P_bar(2,:)';     J2(3:3:(3*numP), 2) = -P_bar(1,:)';
    g(1:2:(2*numP)) = P_hat(1,:)'- targetPC_2D(1,:)';
    g(2:2:(2*numP)) = P_hat(2,:)'- targetPC_2D(2,:)';
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
    
    c2 = -kron(((Var.f*P_bar(2,:)')./((Var.f - P_bar(3,:)').^2))/Var.s,ones(1,6));
    J(2:2:(2*numP),3:8) = J2(2:3:(3*numP),:).*c1 + J2(3:3:(3*numP),:).*c2;
    
    A = J'*W*J;
    b = J'*W*g;
    x = (A(2:8,2:8)+eye(7)*2e-2)\b(2:8);
  
    if norm(x) < 1e-6
        break
    end
    
    energy = obj_val(sourcePC_3D, targetPC_2D, Var);
    % Perform line search
    alpha = 1;
    success = 0;
    for searchIter = 1:10
        Var1.s = Var.s;
        Var1.f = Var.f + x(1)*alpha;
        Var1.o = Var.o + x(5:7)*alpha;
        C = [0, -x(4), x(3);
            x(4), 0, -x(2);
            -x(3), x(2), 0]*alpha;
        Var1.R = normalize_rot(expm(C)*Var.R);
        e = obj_val(sourcePC_3D, targetPC_2D, Var1);
        if e < energy
            Var = Var1;
            success = 1;
        end
        alpha = alpha/2;
    end
    if success == 0
        break;
    end
end

Camera_opt = Camera_init;
Camera_opt.scale = Var.s;
Camera_opt.lookAt = Var.o;
Camera_opt.origin = Var.o + Var.R(3,:)'*Var.f;
Camera_opt.upVec = Var.R(2,:)';
d_ori = norm(Camera_opt.lookAt-Camera_init.lookAt);
d_eye = norm(Camera_opt.origin-Camera_init.origin);
d_upVec = norm(Camera_opt.upVec-Camera_init.upVec);
fprintf('d_ori = %f, d_eye = %f, d_upVec = %f.\n', d_ori, d_eye, d_upVec); 

function [e] = obj_val(sourcePC_3D, targetPC_2D, Var)
%
numP = size(sourcePC_3D, 2);
P_bar = Var.R*(sourcePC_3D - Var.o*ones(1,numP));
P_hat(1,:) = Var.f*P_bar(1,:)./(Var.f - P_bar(3,:));
P_hat(2,:) = Var.f*P_bar(2,:)./(Var.f - P_bar(3,:));
P_hat = P_hat/Var.s;
g(1:2:(2*numP)) = P_hat(1,:)'- targetPC_2D(1,:)';
g(2:2:(2*numP)) = P_hat(2,:)'- targetPC_2D(2,:)';
e = g*g';

function [R] = normalize_rot(R_in)
R_in(:,1) = R_in(:,1)/norm(R_in(:,1));
R_in(:,2) = R_in(:,2) - R_in(:,1)*(R_in(:,2)'*R_in(:,1));
R_in(:,2) = R_in(:,2)/norm(R_in(:,2));
R_in(:,3) = R_in(:,3) - R_in(:,1)*(R_in(:,3)'*R_in(:,1)) - R_in(:,2)*(R_in(:,3)'*R_in(:,2));
R_in(:,3) = R_in(:,3)/norm(R_in(:,3));
R = R_in;
