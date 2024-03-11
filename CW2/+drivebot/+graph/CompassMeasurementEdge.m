classdef CompassMeasurementEdge < g2o.core.BaseUnaryEdge
   
    % Q1c:
    % This implementation contains a bug. Identify the problem
    % and fix it as per the question.

    properties(Access = protected)
        
        compassAngularOffset;
        
    end
    
    methods(Access = public)
    
        function this = CompassMeasurementEdge(compassAngularOffset)
            this = this@g2o.core.BaseUnaryEdge(1);
            this.compassAngularOffset = compassAngularOffset;
        end
        
        function computeError(this)
            x = this.edgeVertices{1}.estimate();
            % wrap the angle between -pi to pi 
            this.errorZ =  g2o.stuff.normalize_theta(x(3) + this.compassAngularOffset - this.z);
        end
        
        function linearizeOplus(this)
            % Compute Jacobian matrix which is the partial derivative of the error function with respect to X,Y,angle  
            this.J{1} = [0 0 1];
        end        
    end
end