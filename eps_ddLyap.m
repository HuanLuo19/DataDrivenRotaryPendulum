classdef eps_ddLyap < handle
    %DDLYAP 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        x ,
        l {mustBeNumeric, mustBeFinite},
        dtau {mustBeNumeric, mustBeFinite},
        Ki {mustBeNumeric, mustBeFinite},
        sys,
        Q {mustBeNumeric, mustBeFinite},
        R {mustBeNumeric, mustBeFinite},
        B {mustBeNumeric, mustBeFinite},
        e {mustBeNumeric, mustBeFinite},
        
        X , % cell of data
        Pi {mustBeNumeric, mustBeFinite},
        Pi_e {mustBeNumeric, mustBeFinite},
        Pim1 {mustBeNumeric, mustBeFinite},
        Kip1 {mustBeNumeric, mustBeFinite},
        Kip1_e {mustBeNumeric, mustBeFinite},
        Pi_lyap {mustBeNumeric, mustBeFinite},
        Kip1_lyap {mustBeNumeric, mustBeFinite},
        Pi_lyap_e {mustBeNumeric, mustBeFinite},
        Kip1_lyap_e {mustBeNumeric, mustBeFinite},

        data_matrix_end_point , % 
        data_array_integral , %
    end
    
    methods
        function self = eps_ddLyap(x, l, dtau, Pim1, Ki, sys, e)
            %DDLYAP 构造此类的实例
            %   此处显示详细说明
            self.x = x;
            self.dtau = dtau;
            % l: number of scaler equations
            % To use continuous data, set l as a none-zero interger
            % To use seperate data, selt l=0
            self.l = l; 
            self.Pim1 = Pim1;
            self.Ki = Ki;
            self.sys = sys;
            self.Q = sys.Q;
            self.R = sys.R;
            self.B = sys.B;
            self.e = e;

            self.Pi_lyap = lyap((sys.A - self.B * self.Ki)', self.Q + self.Ki' * self.R * self.Ki);
            self.Kip1_lyap = self.R \ self.B' * self.Pi_lyap;

            self.Pi_lyap_e = lyap((sys.A - self.B * self.Ki)',...
                ( -(1-self.e)*((sys.A - self.B * self.Ki)' * self.Pim1 + self.Pim1 * (sys.A - self.B * self.Ki)) +...
                self.e*(self.Q + self.Ki' * self.R * self.Ki)));
            self.Kip1_lyap_e = self.R \ self.B' * self.Pi_lyap_e;
            
            if self.l
                self.solveDataEqu_contData();
            else
                self.solveDataEqu_sepData();
            end
        end

        function solveDataEqu_contData(self)
            % When passing one continuos data trajectory, this fcn will
            % evenly divide the data to 'l' segments to formulate scaler
            % equations

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
                delta_xx(i,:) = self.traj2bar(self.X{i}(:,end))' - self.traj2bar(self.X{i}(:,1))'; % delta_xx = [x(t1) - x(t0),...]
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
            % vecs_Pi = lsqminnorm(delta_xx, - I_xx * vec_Qk);   % calculate using least-square
            % vecs_Pi = pinv(delta_xx) * - I_xx * vec_Qk;          % calculate using pseudo-inverse
            vecs_Pi = pinv(delta_xx) * ( (1-self.e) * delta_xx * vecs_Pim1 - self.e * I_xx * vec_Qk);          % calculate using pseudo-inverse
            self.Pi = self.vecs2mat(vecs_Pi);             % ddLyap solution

            self.Kip1 = self.R \ self.B' * self.Pi;     % policy iteration

            % --- test
            self.data_matrix_end_point = delta_xx;
            self.data_array_integral = I_xx;
            % --- end test
        end

        function solveDataEqu_sepData(self)
            % When passing seperate data trajectory as cell, this fcn will
            % use each row of data in the cell to formulate scaler equations

            % define data segments
            self.X = self.x;
            N = size(self.X,1);
            n = size(self.X{1},1);
            
            % compute delta_xx
            delta_xx = zeros(N, n*(n+1)/2);
            for i = 1 : N
                delta_xx(i,:) = self.traj2bar(self.X{i}(:,end))' - self.traj2bar(self.X{i}(:,1))';
            end
            
            % compute Ixx
            XkX = cell(N,1);
            int_XkX = cell(N,1);
            I_xx = zeros(N, n*n);
            for j = 1:N
                XkX{j} = zeros(n*n, size(self.X{j},2));
                for i = 1:size(self.X{1},2)
                    XkX{j}(:,i) = kron(self.X{j}(:,i),self.X{j}(:,i)); % xx(t) = x(t) (kron) x(t)
                end
                int_XkX{j} = cumtrapz(self.dtau, XkX{j}, 2);
                I_xx(j,:) = (int_XkX{j}(:,end) - int_XkX{j}(:,1))';
            end
            
            Qk = self.Q + self.Ki' * self.R * self.Ki;
            vec_Qk = reshape(Qk,[n*n,1]);

            % step size solution
            vecs_Pi_e = pinv(delta_xx) * ( (1-self.e) * delta_xx * self.mat2vecs(self.Pim1) - self.e * I_xx * vec_Qk);          % calculate using pseudo-inverse
            self.Pi_e = self.vecs2mat(vecs_Pi_e);             % ddLyap solution
            self.Kip1_e = self.R \ self.B' * self.Pi_e;     % policy iteration

            % original solution
            vecs_Pi = pinv(delta_xx) * ( - I_xx * vec_Qk);          % calculate using pseudo-inverse
            self.Pi = self.vecs2mat(vecs_Pi);             % ddLyap solution
            self.Kip1 = self.R \ self.B' * self.Pi;     % policy iteration

            % --- test
            self.data_matrix_end_point = delta_xx;
            self.data_array_integral = I_xx;
            % --- end test
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
            % VECS2MAT
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

        function [vecs_P] = mat2vecs(~,P)
            % MAT2VECS
            % P must be symmetric
            % P = [a b;
            %      b c]
            % vecs_P = [a, 2b, c]'
            nx = size(P,2);
            n = nx * (nx + 1) / 2;
            vecs_P = zeros(n,1);
            temp = 1;
            for i = 1:nx
                for j = i:nx
                    if i==j
                        vecs_P(temp) = P(i,j);
                    else
                        vecs_P(temp) = 2 * P(i,j);
                    end
                    temp = temp + 1;
                end
            end
        end
    
    end
end

