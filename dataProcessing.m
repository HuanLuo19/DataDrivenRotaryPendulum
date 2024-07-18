classdef dataProcessing < handle
    %DATAPROCESSING 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        x , % pre-processed data
        t , % time label
    end
    
    methods
        function self = dataProcessing(x,t)
            %DATAPROCESSING 构造此类的实例
            %   此处显示详细说明
            self.x = x;
            self.t = t;
        end
        
        function [x_out, t_out] = getTimeIntervalData(self,T)
            idx_min = find(abs(self.t - T(1)) < 0.001);
            idx_max = find(abs(self.t - T(2)) < 0.001);
            if isempty(idx_min) || isempty(idx_max)
                fprintf("Cannot find data on the interval end point! \n")
                return
            end
            x_out = self.x(:,idx_min:idx_max);
            t_out = self.t(:,idx_min:idx_max);
        end

        function X = dataDivider(self)
            % this method divide data into desired segments
            % X is the output data in cell form
            % every row of cell X represents a data segments
            % which will be used in the data equation as one scaler
            % equation
        end

        function dataLoss(self)
            % to be continued
        end
    end
end

