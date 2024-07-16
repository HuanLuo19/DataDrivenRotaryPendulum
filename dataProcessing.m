classdef dataProcessing
    %DATAPROCESSING 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        x , % pre-processed data
        Property1
    end
    
    methods
        function self = dataProcessing(x,inputArg2)
            %DATAPROCESSING 构造此类的实例
            %   此处显示详细说明
            self.x = x;
            self.Property1 = x + inputArg2;
        end
        
        function X = dataDivider(obj,inputArg)
            % this method divide data into desired segments
            % X is the output data in cell form
            % every row of cell X represents a data segments
            % which will be used in the data equation as one scaler
            % equation
            X = obj.Property1 + inputArg;
        end
    end
end

