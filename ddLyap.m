classdef ddLyap < handle
    %DDLYAP 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        x {mustBeNumeric, mustBeFinite},
        dtau {mustBeNumeric, mustBeFinite},
        l {mustBeNumeric, mustBeFinite},
        Ki {mustBeNumeric, mustBeFinite},
        sys,
        Q {mustBeNumeric, mustBeFinite},
        R {mustBeNumeric, mustBeFinite},
        B {mustBeNumeric, mustBeFinite},
        
        X , % cell of data
        Pi {mustBeNumeric, mustBeFinite},
        Kip1 {mustBeNumeric, mustBeFinite},
        Pi_lyap {mustBeNumeric, mustBeFinite},
        Kip1_lyap {mustBeNumeric, mustBeFinite},
    end
    
    methods
        function self = ddLyap(x, dtau, l, Ki, sys)
            %DDLYAP 构造此类的实例
            %   此处显示详细说明
            self.x = x;
            self.dtau = dtau;
            self.l = l;
            self.Ki = Ki;
            self.sys = sys;
            self.Q = sys.Q;
            self.R = sys.R;
            self.B = sys.B;

            self.Pi_lyap = lyap((sys.A - self.B * self.Ki)', self.Q + self.Ki' * self.R * self.Ki);
            self.Kip1_lyap = self.R \ self.B' * self.Pi_lyap;
            self.solveDataEqu();
            % self.solveDataEqu_test();
        end

        function solveDataEqu(self)
            %
            % define data segments
            index_int = (size(self.x,2) - 1) / self.l;
            self.X = cell(self.l,1);
            for j = 1:self.l 
                self.X{j} = self.x(:, (index_int*(j-1) + 1 ): (index_int*j + 1 ));
            end

            % compute data matrices
            n = size(self.x,1);
            % compute delta_xx
            delta_xx = zeros(self.l, n*(n+1)/2);
            for i = 1 : self.l
                delta_xx(i,:) = self.traj2bar(self.X{i}(:,end))' - self.traj2bar(self.X{i}(:,1))';
            end
            % compute Ixx
            XkX = cell(self.l,1);
            int_XkX = cell(self.l,1);
            I_xx = zeros(self.l, n*n);
            for j = 1:self.l
                XkX{j} = zeros(n*n, size(self.X{j},2));
                for i = 1:size(self.X{1},2)
                    XkX{j}(:,i) = kron(self.X{j}(:,i),self.X{j}(:,i)); % xx(t) = x(t) (kron) x(t)
                end
                int_XkX{j} = cumtrapz(self.dtau, XkX{j}, 2);
                I_xx(j,:) = (int_XkX{j}(:,end) - int_XkX{j}(:,1))';
            end
            
            Qk = self.Q + self.Ki' * self.R * self.Ki;
            vec_Qk = reshape(Qk,[n*n,1]);
            % Pi_hat = lsqminnorm(delta_xx, - I_xx * vec_Qk);   % calculate using least-square
            vecs_Pi = pinv(delta_xx) * - I_xx * vec_Qk;          % calculate using pseudo-inverse
            self.Pi = self.vecs2mat(vecs_Pi);             % ddLyap solution

            self.Kip1 = self.R \ self.B' * self.Pi;     % policy iteration
        end

            function solveDataEqu_test(self)
            % test20240628
            % define data segments
            index_int = (size(self.x,2) - 1) / self.l;
            self.X = cell(self.l,1);
            for j = 1:self.l 
                self.X{j} = self.x(:, (index_int*(j-1) + 1 ): (index_int*j + 1 ));
            end

            % compute data matrices
            n = size(self.x,1);
            % compute delta_xx
            delta_xx = zeros(self.l, n*n);
            for i = 1 : self.l
                delta_xx(i,:) = kron(self.X{i}(:,end), self.X{i}(:,end))' ...
                                - kron(self.X{i}(:,1), self.X{i}(:,1))';
            end

            % compute Ixx
            XkX = cell(self.l,1);
            int_XkX = cell(self.l,1);
            I_xx = zeros(self.l, n*n);
            for j = 1:self.l
                XkX{j} = zeros(n*n, size(self.X{j},2));
                for i = 1:size(self.X{1},2)
                    XkX{j}(:,i) = kron(self.X{j}(:,i),self.X{j}(:,i)); % xx(t) = x(t) (kron) x(t)
                end
                int_XkX{j} = cumtrapz(self.dtau, XkX{j}, 2);
                I_xx(j,:) = (int_XkX{j}(:,end) - int_XkX{j}(:,1))';
            end
            
            Q_k = self.Q + self.Ki' * self.R * self.Ki;
            varXi_k = - I_xx * reshape(Q_k,[n*n, 1]);
            % Pi_hat = lsqminnorm(delta_xx, varXi_k);   % calculate using least-square
            vecs_Pi = pinv(delta_xx) * varXi_k;          % calculate using pseudo-inverse
            self.Pi = reshape(vecs_Pi,[n, n]);             % ddLyap solution

            self.Kip1 = self.R \ self.B' * self.Pi;     % policy iteration
        end

        function [x_bar] = traj2bar(~,x)
            % traj2bar 
            % x = [x1, x2]'
            % x_bar = [x1^2, x1x2, x2^2]';
            nx = size(x,1);
            N = nx * (nx + 1) / 2;
            
            x_bar = zeros(N,size(x,2));
            for col = 1:size(x,2)
                row = 1;
                for i = 1:nx
                    for j = i:nx
                        x_bar(row,col) = x(i,col) * x(j,col);
                        row = row + 1;
                    end
                end
            end
        end

        function [P] = vecs2mat(~,vecs_P)
            % HAT2MAT
            % P must be symmetric
            % P = [a b;
            %      b c]
            % vecs_P = [a, 2b, c]'
            n = (-1 + sqrt(1 + 8*size(vecs_P,1))) / 2;
            P = zeros(n,n);
            temp = 1;
            for i = 1:n
                for j = i:n
                    if i==j
                        P(i,j) = vecs_P(temp);
                    else
                        P(i,j) = vecs_P(temp) / 2;
                        P(j,i) = P(i,j);
                    end
                    temp = temp + 1;
                end
            end
        end
    
    end
end

