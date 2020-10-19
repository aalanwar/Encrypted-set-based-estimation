function updatePlotConZono(varargin)

%If only two  argument is passed
if nargin==2
    h = varargin{1};
    Z=varargin{2};
    dims=[1,2];
    type{1}='b';
    
%If three arguments are passed    
elseif nargin==3
    h = varargin{1};
    Z=varargin{2};
    dims=varargin{3};
    type{1}='b';
    
%If four or more arguments are passed
elseif nargin>=4
    h = varargin{1};
    Z=varargin{2};
    dims=varargin{3};   
    type(1:length(varargin)-3)=varargin(4:end);
end

% project zonotope
Z = project(Z,dims);

% delete zero generators
if isempty(Z.A)
    p = polygon(zonotope(Z.Z));
    %plot and output the handle
    set(h,'XData', p(1,:), 'YData', zeros(size(p(1,:))),'ZData',p(2,:));
else
    p = vertices(Z);
    %add the first vertice to the end to close the shape
    ver = p;
    ver(:,end+1)=ver(:,1);
    %plot and output the handle
    set(h,'XData', ver(1,:), 'YData', zeros(size(ver(1,:))),'ZData',ver(2,:));
end



%------------- END OF CODE --------------