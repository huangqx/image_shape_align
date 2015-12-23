function [R, t] = horn87(P, Q, w)

P = double(P);
Q = double(Q);
w = double(w);

nP = size(P, 2);

W = sparse(1:length(w), 1:length(w), w);

P1 = P*W;
Q1 = Q*W;

cP = sum(P1')'/sum(w);
cQ = sum(Q1')'/sum(w);

P = P - cP*ones(1, nP);
Q = Q - cQ*ones(1, nP);

S = P*W*Q';

N = zeros(4,4);
N(1,1) = S(1,1) + S(2,2) + S(3,3);
N(1,2) = S(2,3) - S(3,2);
N(1,3) = S(3,1) - S(1,3);
N(1,4) = S(1,2) - S(2,1);
N(2,1) = N(1,2);
N(3,1) = N(1,3);
N(4,1) = N(1,4);

N(2,2) = S(1,1) - S(2,2) - S(3,3);
N(2,3) = S(1,2) + S(2,1);
N(2,4) = S(3,1) + S(1,3);
N(3,2) = N(2,3);
N(4,2) = N(2,4);

N(3,3) = S(2,2) - S(1,1) - S(3,3);
N(3,4) = S(2,3) + S(3,2);
N(4,3) = N(3,4);

N(4,4) = S(3,3) - S(1,1) - S(2,2);

[u,v] = eig(N);

q = u(:,4);
q0 = q(1);
n = q(2:4);
cor = [0, -n(3), n(2);
    n(3), 0, -n(1);
    -n(2), n(1), 0];

R = (q0*q0-n'*n)*eye(3) + 2*n*n' + 2*q0*cor;
t = cQ - R*cP;
