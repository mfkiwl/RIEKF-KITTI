clear;clc;

mat = randn(3, 3);
[u, s, v] = svd(mat);

R = u;
% SO(3) -> so(3)
theta = acos((trace(R)-1)/2);
lnR = theta/(2*sin(theta))*(R - R');
t1 = skew2vec(lnR);

g0 = Gamma0(t1);
g1 = Gamma1(t1);
g2 = Gamma2(t1);

g0_ = GammaFunc(0, t1);
g1_ = GammaFunc(1, t1);
g2_ = GammaFunc(2, t1);

g = [0, 0, -9.81];
A_r = zeros(9, 9);
A_r(4:6, 1:3) = vec2skew(g);
A_r(7:9, 4:6) = eye(3);

filter = RI_EKF();
filter.prediction(randn(6, 1));
filter.propagation([1;2;3]);

function vec = skew2vec(mat)
    vec = zeros(3, 1);
    vec(1) = mat(3, 2);
    vec(2) = mat(1, 3);
    vec(3) = mat(2, 1);
end

function [Ax] = vec2skew(v)
% Convert from vector to skew symmetric matrix
Ax = [    0, -v(3),  v(2);
       v(3),     0, -v(1);
      -v(2),  v(1),     0];
end

function val = Gamma0(twist)
% exponential map from so(3) to SO(3).
    theta = sqrt(twist'*twist);
    twist_ = vec2skew(twist);
    val = eye(3) + (sin(theta)/theta)*twist_ + ...
          ((1-cos(theta))/theta^2)*twist_^2;
end

function val = Gamma1(twist)
% left Jacobian of SO(3).
    theta = norm(twist, 2);
    twist_ = vec2skew(twist);
    val = eye(3) + ((1-cos(theta))/theta^2) * twist_ + ...
          ((theta - sin(theta))/theta^3)*twist_^2;
end

function val = Gamma2(twist)
    theta = norm(twist, 2);
    twist_ = vec2skew(twist);
    val = 0.5*eye(3) + ((theta-sin(theta))/theta^3)*twist_...
        + ((theta^2 + 2*cos(theta) - 2)/(2*theta^4))*twist_^2;
end

function val = GammaFunc(m, twist)
    num = 10;
    val = zeros(3, 3);
    twist_ = vec2skew(twist);
    for n = 0:num
       val = val + twist_^n/factorial(n+m);
    end
end