 classdef RI_EKF < handle
    properties
        mu;          % state;
        Sigma;       % covariance;
        delta_t;     % time step;
        g;           % gravity vector;
        vk;          % measurement noise;
        Q;           % IMU noise;
        bias;        % IMU bias;
    end
   
    methods
        function obj = RI_EKF(mu_init, sigma_init, v_init)
            % state:
            % X = [R, v, p;
            %      0, 1, 0;
            %      0, 0, 1]
            obj.mu = eye(5);
            obj.mu(1:3, 1:3) = mu_init(1:3, 1:3);
            obj.mu(1:3, 5) = mu_init(1:3, 4);
            obj.mu(1:3, 4) = v_init;
            obj.Sigma = sigma_init;
            
            obj.delta_t = 0.1;
            obj.g = [0; 0; -9.81];
            obj.vk = blkdiag(0.1*eye(3), zeros(2));
            obj.Q  = 0.1*eye(15); obj.Q(7:9, 7:9) = zeros(3);
            obj.bias = zeros(6, 1);
        end
        
        function prediction(obj, u)
            phi = obj.Phi();
            % input: 6-vector: [omega; a].
            u = u - obj.bias;
            w = u(1:3);
            a = u(4:6);
            w = obj.vec2skew(w);
            
            % Update state with IMU model.
            R = obj.mu(1:3, 1:3);
            v = obj.mu(1:3, 4);
            p = obj.mu(1:3, 5);
            AdX = obj.Adj(obj.mu);
            
            obj.mu(1:3, 1:3) = R*obj.Gamma(0, w*obj.delta_t);
            obj.mu(1:3, 4) = v + R*obj.Gamma(1, w*obj.delta_t)*a*obj.delta_t + obj.g*obj.delta_t;
            obj.mu(1:3, 5) = p + v*obj.delta_t + R*obj.Gamma(2, w*obj.delta_t)*a*obj.delta_t^2 + 0.5*obj.g*obj.delta_t^2;        
        
            % Covariance propagation.
            AdX = blkdiag(AdX, eye(6));
            obj.Sigma = phi*obj.Sigma*phi' + obj.delta_t*AdX*phi*obj.Q*phi'*AdX';
        end
        
        function correction(obj, y)
            % input: y = [x, y, z], GPS information;
            y = [y; 0; 1];
            H = [zeros(3), zeros(3), -eye(3), zeros(3), zeros(3)];
            b = [zeros(4, 1); 1];
            R = obj.mu(1:3, 1:3);
            N = R * obj.vk(1:3, 1:3) * R';
            
            S = H*obj.Sigma*H' + N(1:3, 1:3);
%             disp(det(S))
            K = obj.Sigma*H' / S;
            K_s = K(1:9, :);
            K_b = K(10:15, :);
            
            PI = [eye(3), zeros(3, 2)];
            nv = obj.mu*y - b;
            
            obj.bias = obj.bias + K_b*PI*nv;
            obj.mu = expm(obj.stateMat(K_s*PI*nv))*obj.mu;
            obj.Sigma = (eye(15)-K*H)*obj.Sigma*(eye(15)-K*H)' + K*N(1:3, 1:3)*K';
        end
        
        function phi = Phi(obj)
            R = obj.mu(1:3, 1:3);
            v = obj.mu(1:3, 4);
            p = obj.mu(1:3, 5);
            
            A_r = zeros(15, 15);
            A_r(4:6, 1:3) = obj.vec2skew(obj.g);
            A_r(7:9, 4:6) = eye(3);
            A_r(1:3, 10:12) = -R;
            A_r(4:6, 10:12) = -obj.vec2skew(v)*R;
            A_r(4:6, 13:15) = -R;
            A_r(7:9, 10:12) = -obj.vec2skew(p)*R;

            phi = expm(A_r*obj.delta_t);            
        end
        
        function state = stateMat(obj, vec)
            % input: vec = [9x1].
            % including rotation, velocity, position, IMU bias.
            twistR = vec(1:3);
            v = vec(4:6);
            p = vec(7:9);
            state = blkdiag(eye(3), zeros(2));
            state(1:3, 1:3) = obj.vec2skew(twistR);
            state(1:3, 4) = v;
            state(1:3, 5) = p;
        end

        function R = euler2so3(obj, w)
            % input: w = [roll, pitch, yaw]
            rx = w(1);
            ry = w(2);
            rz = w(3);
            Rx = [1 0 0; 0 cos(rx) -sin(rx); 0 sin(rx) cos(rx)]; % base => nav  (level oxts => rotated oxts)
            Ry = [cos(ry) 0 sin(ry); 0 1 0; -sin(ry) 0 cos(ry)]; % base => nav  (level oxts => rotated oxts)
            Rz = [cos(rz) -sin(rz) 0; sin(rz) cos(rz) 0; 0 0 1];
            R  = Rz*Ry*Rx;
            R = obj.skew2vec(logm(R));
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