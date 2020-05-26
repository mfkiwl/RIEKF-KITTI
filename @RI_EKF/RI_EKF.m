classdef RI_EKF < handle
    properties
        mu;          % state;
        Sigma;       % covariance;
        gfun;        % dynamic function;
        zfun;        % observation function;
        delta_t;     % time step;
        g;           % gravity vector;
        vk;          % measurement noise;
        Q;           % IMU noise;
    end
   
    methods
        function obj = RI_EKF()
            % state:
            % X = [R, v, p;
            %      0, 1, 0;
            %      0, 0, 1]
            obj.mu = eye(5);
            obj.Sigma = 0.1*eye(9);
            obj.delta_t = 1;
            obj.g = [0; 0; -9.81];
            obj.vk = eye(5);
            obj.Q  = 100*eye(9);
        end
        
        function prediction(obj, u)
            % input: 6-vector: [omega; a].
            omega = u(1:3);
            a = u(4:6);
            w = obj.euler2SO3(omega);
            
            % Update state with IMU model.
            R = obj.mu(1:3, 1:3);
            v = obj.mu(1:3, 4);
            p = obj.mu(1:3, 5);
            AdX = obj.Adj(obj.mu);
            obj.mu(1:3, 1:3) = R*obj.Gamma(0, logm(w)*obj.delta_t);
            obj.mu(1:3, 4) = v + R*obj.Gamma(1, logm(w)*obj.delta_t)*a + obj.g*obj.delta_t;
            obj.mu(1:3, 5) = p + v*obj.delta_t + R*obj.Gamma(2, logm(w)*obj.delta_t)*a*obj.delta_t^2 + 0.5*obj.g*obj.delta_t^2;        
        
            % Covariance propagation.
            A_r = zeros(9, 9);
            A_r(4:6, 1:3) = obj.vec2skew(obj.g);
            A_r(7:9, 4:6) = eye(3);
            % Phi: transition matrix.
            Phi = expm(A_r);
            obj.Sigma = Phi*obj.Sigma*Phi' + AdX*obj.Q*AdX';
        end
        
        function propagation(obj, y)
            % input: y = [x, y, z], GPS information;
            y = [y; 0; 1];
            H = [zeros(3, 3), zeros(3, 3), eye(3)];
            b = [zeros(4, 1); 1];
            N = obj.mu * obj.vk * obj.mu';
            S = H*obj.Sigma*H' + N(1:3, 1:3);
            L = obj.Sigma*H' / S;
            
            r = obj.mu*y - b;
            r = r(1:3);
            r = L * r;
            obj.mu = expm(obj.stateMat(r))*obj.mu;
            obj.Sigma = (eye(9)-L*H)*obj.Sigma*(eye(9)-L*H)' + L*N(1:3, 1:3)*L';
        end
        
        function state = stateMat(obj, vec)
            % input: vec = [9x1].
            twistR = vec(1:3);
            v = vec(4:6);
            p = vec(7:9);
            state = eye(5);
            state(1:3, 1:3) = obj.vec2skew(twistR);
            state(1:3, 4) = v;
            state(1:3, 5) = p;
        end
        
        function R = euler2SO3(obj, w)
            % input: w = [roll, pitch, yaw]
            rx = w(1);
            ry = w(2);
            rz = w(3);
            Rx = [1 0 0; 0 cos(rx) -sin(rx); 0 sin(rx) cos(rx)]; % base => nav  (level oxts => rotated oxts)
            Ry = [cos(ry) 0 sin(ry); 0 1 0; -sin(ry) 0 cos(ry)]; % base => nav  (level oxts => rotated oxts)
            Rz = [cos(rz) -sin(rz) 0; sin(rz) cos(rz) 0; 0 0 1];
            R  = Rz*Ry*Rx;
        end
       
        function val_g = Gamma(obj, m, twist)
            num = 10;
            val_g = zeros(3, 3);
            twist_ = obj.vec2skew(twist);
            for n = 0:num
                val_g = val_g + twist_^n/factorial(n+m);
            end
        end
        
        function AdX = Adj(obj, X)
            % adjoint for SE_2(3).
            R = X(1:3, 1:3);
            v = X(1:3, 4);
            p = X(1:3, 5);
            AdX = blkdiag(R, R, R);
            AdX(4:6, 1:3) = obj.vec2skew(v)*R;
            AdX(7:9, 1:3) = obj.vec2skew(p)*R;
        end
       
        function skew = vec2skew(obj, vec)
            skew = [    0, -vec(3),  vec(2);
                   vec(3),       0, -vec(1);
                  -vec(2),  vec(1),      0];
        end
       
        function vec = skew2vec(obj, skew)
            vec = zeros(3, 1);
            vec(1) = skew(3, 2);
            vec(2) = skew(1, 3);
            vec(3) = skew(2, 1);    
        end
   end
end