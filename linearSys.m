classdef linearSys < handle
    % Class of LQR problem
    % system: \dot x = Ax + Bu
    % cost matrices: Q,R
    
    properties (Access = public)
        % defined by user
        A {mustBeNumeric, mustBeFinite};
        B {mustBeNumeric, mustBeFinite};
        x0 {mustBeNumeric, mustBeFinite};
        Q {mustBeNumeric, mustBeFinite};
        R {mustBeNumeric, mustBeFinite};

        % calculated from system
        n;          % dimension of state
        m;          % dimension of input
        K_opt;      % optimal gain K_opt calculated from lqr()
        P_opt;      % optimal P_opt calculated from lqr()
    end

    properties (Access = private)
        dtau_sample = 1e-4;
    end
    
    methods
        function self = linearSys(A,B,x0,Q,R)
            % Initializing object
            self.A = A;
            self.B = B;
            self.x0 = x0;
            self.n = size(A,2);
            self.m = size(B,2);

            self.Q = Q;
            self.R = R;
            [self.K_opt, self.P_opt, ~] = lqr(A,B,Q,R,[]);
        end

        function [x,u,tspan] = ResponseFromGain(self,K,dtau,T)
            % generate data in [0,T] with u(t) = -K*x(t)
            if size(K,1) ~= self.m
                fprintf('rows of K not correct \n');
            elseif size(K,2) ~= self.n
                fprintf('columns of K not correct \n');
            end
            tspan = 0:dtau:T;
            x_cell = arrayfun(@(t) expm((self.A - self.B * K)*t)*self.x0, tspan, ...
                'UniformOutput', false);
            x = cell2mat(x_cell);

            u = K * x;

            % % plot
            % figure
            % % sgtitle("K = " + num2str(K),'Interpreter','latex')
            % % subplot(2,1,1)
            % plot(tspan,x,'.',MarkerSize=10);
            % xlabel('$t$','Interpreter','latex')
            % ylabel('$x$','Interpreter','latex')
            % % ylim([-0.3 0.3])
            % ax = gca;
            % ax.FontSize = 14;
            % % 
            % % subplot(2,1,2)
            % % plot(tspan,u,'.',MarkerSize=10);
            % % xlabel('$t$','Interpreter','latex')
            % % ylabel("$u=Kx$",'Interpreter','latex')
            % % ax = gca;
            % % ax.FontSize = 14;
        end

        function generateTestResponseZeroInput(self,T)
            % generate data in [0,T] with u(t) = 0
            tspan = 0:self.dtau_sample:T;
            x_free_cell = arrayfun(@(t) expm(self.A*t)*self.x0, tspan, ...
                'UniformOutput', false);
            x_free = cell2mat(x_free_cell);

            % plot
            figure
            plot(tspan,x_free,LineWidth=2);
            xlabel('$t$','Interpreter','latex')
            title('$x$','Interpreter','latex')
            ax = gca;
            ax.FontSize = 14;
        end

        function generateTestResponseOptimalInput(self,T)
            % generate data in [0,T] with u(t) = -K_opt*x(t)
            tspan = 0:self.dtau_sample:T;
            x_opt_cell = arrayfun(@(t) expm((self.A - self.B * self.K_opt)*t)*self.x0, tspan, ...
                'UniformOutput', false);
            x_opt = cell2mat(x_opt_cell);

            u_opt = self.K_opt * x_opt;

            % plot
            figure
            subplot(2,1,1)
            plot(tspan,x_opt,'.',MarkerSize=10);
            xlabel('$t$','Interpreter','latex')
            title('$x$','Interpreter','latex')
            ax = gca;
            ax.FontSize = 14;
            
            subplot(2,1,2)
            plot(tspan,u_opt,'.',MarkerSize=10);
            xlabel('$t$','Interpreter','latex')
            title('$u^*$','Interpreter','latex')
            ax = gca;
            ax.FontSize = 14;
        end
    end
end

