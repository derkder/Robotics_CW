% This edge encodes a 3D range bearing measurement.
%
% The measurement is in spherical polar coordinates

% Jacobian from https://github.com/petercorke/robotics-toolbox-matlab/blob/master/RangeBearingSensor.m

% 这里的公式和workshop2里的公式一摸一样,处理beita_k的最后一项
classdef LandmarkRangeBearingEdge < g2o.core.BaseBinaryEdge
    properties(Access = public)
    end


    methods(Access = public)
    
        function this = LandmarkRangeBearingEdge()
            this = this@g2o.core.BaseBinaryEdge(2);
        end
        
        function initialize(this)
            % Q2b:
            % Complete implementation
            x = this.edgeVertices{1}.x;
            deltaX = this.z(1) * cos(this.z(2) + x(3));
            deltaY = this.z(1) * sin(this.z(2) + x(3));
            landmark = [x(1) + deltaX; x(2) + deltaY];
            % disp(landmark);
            this.edgeVertices{2}.setEstimate(landmark);
            
        end
        
        % [this.range, this.bearing]
        function computeError(this)
            x = this.edgeVertices{1}.x;
            %disp(x)
            %disp(this.edgeVertices{2}.estimate())
            dx = this.edgeVertices{2}.estimate() - this.edgeVertices{1}.x([1 2]);
            
            
            this.errorZ(1) = norm(dx) - this.z(1);
            this.errorZ(2) = g2o.stuff.normalize_theta(atan2(dx(2), dx(1)) - x(3) - this.z(2));
        end
        
        function linearizeOplus(this) 
            %x = this.edgeVertices{1}.estimate();
            dx = this.edgeVertices{2}.x - this.edgeVertices{1}.x([1 2]);
            r = norm(dx);
            
            this.J{1} = ...
                [-dx(1)/r -dx(2)/r 0;
                dx(2)/r^2 -dx(1)/r^2 -1];
            this.J{2} = - this.J{1}(1:2, 1:2);
        end        
    end
end