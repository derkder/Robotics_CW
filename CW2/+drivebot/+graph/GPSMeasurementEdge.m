classdef GPSMeasurementEdge < g2o.core.BaseUnaryEdge
   
    properties(Access = protected)
        
        xyOffset;
        
    end
    
    methods(Access = public)
    
        function this = GPSMeasurementEdge(xyOffset)
            this = this@g2o.core.BaseUnaryEdge(2);
            this.xyOffset = xyOffset;
        end
        
        function computeError(this)

	        % Q1d:
            % Implement the code
            %warning('gpsmeasurementedge:computeerror:unimplemented', ...
             %       'Implement the rest of this method for Q1d.');
            x = this.edgeVertices{1}.estimate();
            %disp([this.xyOffset, x]);
            c = cos(x(3));
            s = sin(x(3));
            M = [c -s; s c];
            % the error in the GPS observation function which only considers X and Y  
            this.errorZ = (x(1:2) + M * this.xyOffset) - this.z;
        end
        
        function linearizeOplus(this)

	    % Q1d:
        % Implement the code
        %warning('gpsmeasurementedge:lineareizeoplus:unimplemented', ...
        %        'Implement the rest of this method for Q1d.');
        x = this.edgeVertices{1}.estimate();
        % Compute Jacobian matrix which is the partial derivative of the error function with respect to X,Y,angle  
        c = cos(x(3));
        s = sin(x(3));
        dX = this.xyOffset(1);
        dY = this.xyOffset(2);
        this.J{1} = [-1 0 dX * s + dY * c; 0 -1 -dX * c + dY * s];
        end
    end
end
