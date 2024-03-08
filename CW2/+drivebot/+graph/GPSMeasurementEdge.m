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
        this.errorZ(1) = x(1) + this.xyOffset(1) - this.z(1);
        this.errorZ(2) = x(2) + this.xyOffset(2) - this.z(2);
        end
        
        function linearizeOplus(this)

	    % Q1d:
        % Implement the code
        %warning('gpsmeasurementedge:lineareizeoplus:unimplemented', ...
        %        'Implement the rest of this method for Q1d.');
        x = this.edgeVertices{1}.estimate();

        c = cos(x(3));
        s = sin(x(3));
        dX = this.xyOffset(1);
        dY = this.xyOffset(2);
        this.J{1} = [1 0 -dX * s - dY * c; 0 1 dX * c - dY * s];
        end
    end
end
