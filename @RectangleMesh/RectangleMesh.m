classdef RectangleMesh
    
    properties(SetAccess = private)
        coordinates;
        elements;
    end
    
    % Class constructor
    methods (Access = public, Static = false)
        function self = RectangleMesh(x0, y0, x1, y1, nx, ny)
            % Triangular mesh of the 2D rectangle (x0, y0) x (x1, y1).
            % Given the number of cells (nx, ny) in each direction,
            % the total number of triangles will be 2*nx*ny and
            % the total number of vertices will be (nx + 1)*(ny + 1).
            
            x = linspace(x0, x1, nx + 1);
            y = linspace(y0, y1, ny + 1);
            [xx, yy] = ndgrid(x,y);
            self.coordinates = [xx(:), yy(:)];
            
            %
            elem = [1, 2, nx + 3, nx + 2];
            elem = kron(ones(nx,1), elem) + kron((0:nx-1)', ones(size(elem)));
            elem = kron(ones(ny,1), elem) + kron((0:ny-1)'*(nx+1), ones(size(elem)));
            %
            self.elements = elem;
            
%             fprintf('-- mesh: nb vertices = %d, nb triangles = %d\n', ...
%                 self.getNumCoordinates(), self.getNumElements());
        end
        
        function showMesh(self)
            nc = self.getNumCoordinates();
            trimesh(self.elements, self.coordinates(:, 1), ...
                self.coordinates(:, 2), zeros(nc, 1), 'edgecolor','k'); hold on
            plot(self.coordinates(:,1), self.coordinates(:,2),'ro');
%             plot(center(:,1),center(:,2),'o','LineWidth',1,'MarkerEdgeColor','k','MarkerSize',20);
%             text(center(:,1)-0.025,center(:,2),int2str(elem_range),'FontSize',12,'FontWeight','bold','Color','k');
            text(self.coordinates((1:nc)',1)+0.0125, ...
                 self.coordinates((1:nc)',2)+0.0125, ...
                 int2str((1:nc)'), 'FontSize',10,'FontWeight','bold'); hold off
        
            axis square off;
            view(2);
        end
        
        function output = getNumCoordinates(self)
            output = size(self.coordinates,1);
        end
        
        function output = getNumElements(self)
            output = size(self.elements,1);
        end
    end
end