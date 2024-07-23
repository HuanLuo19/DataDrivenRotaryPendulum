classdef dataProcessing < handle
    %DATAPROCESSING 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        x , % pre-processed data, must be coloum time sequence
        t , % time label, must be coloum time sequence
        switch_state, % switch state data, must be coloum time sequence
    end
    
    methods

        function self = dataProcessing(x,t,switch_state)
            %DATAPROCESSING 构造此类的实例
            %   此处显示详细说明
            self.x = x;
            self.t = t;
            self.switch_state = switch_state;
        end
        
        function [x_on, time_on, switch_state_on] = getSwitchOnData(self)
            idx_switchON_all = find(self.switch_state); % get all index when the switch is on

            time_on = self.t(idx_switchON_all);
            time_on = time_on - time_on(1); % RESET time label from 0
            switch_state_on = self.switch_state(idx_switchON_all);
            x_on = self.x(:,idx_switchON_all);
        end

        function [x_int, t_int, swt_int] = getTimeIntervalData(self,T_int)
            idx_min = find(abs(self.t - T_int(1)) < 0.001);
            idx_max = find(abs(self.t - T_int(2)) < 0.001);
            if isempty(idx_min) || isempty(idx_max)
                fprintf("Cannot find data on the interval end point! \n")
                return
            end
            x_int = self.x(:,idx_min:idx_max);
            t_int = self.t(:,idx_min:idx_max);
            if ~isempty(self.switch_state)
                swt_int = self.switch_state(:,idx_min:idx_max);
            else
                swt_int = [];
            end
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

        function dataSmoothing(self)
            % to be continued
        end
    end
end

