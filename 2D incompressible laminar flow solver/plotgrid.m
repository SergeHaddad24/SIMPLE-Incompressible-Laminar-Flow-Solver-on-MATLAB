function plotgrid(X,Y,cell_centroidX,cell_centroidY,f_centroid_east_x,f_centroid_east_y,f_centroid_west_x,f_centroid_west_y,f_centroid_top_x,f_centroid_top_y,f_centroid_bottom_x,f_centroid_bottom_y)

% plotgrid: To plot structured grid.
%
%    plotgrid(X,Y)
%
%    INPUT:
%      X (matrix)    - matrix with x-coordinates of gridpoints
%      Y (matrix)    - matrix with y-coordinates of gridpoints


if any(size(X)~=size(Y))
   error('Dimensions of X and Y must be equal');
end

[m,n]=size(X);

% Plot grid
figure
set(gcf,'color','w') ;
axis equal
axis off
box on
hold on

% Plot internal grid lines
for i=1:m
    plot(X(i,:),Y(i,:),'b','linewidth',1); 
end
for j=1:n
    plot(X(:,j),Y(:,j),'b','linewidth',1); 
end
plot(cell_centroidX(:),cell_centroidY(:), 'ro','MarkerSize', 2,'MarkerFace','r');
plot(f_centroid_east_x(:),f_centroid_east_y(:),'ro','MarkerSize',3,'MarkerFace','b');
plot(f_centroid_west_x(:),f_centroid_west_y(:),'ro','MarkerSize',3,'MarkerFace','b');
plot(f_centroid_top_x(:),f_centroid_top_y(:),'ro','MarkerSize',3,'MarkerFace','b');
plot(f_centroid_bottom_x(:),f_centroid_bottom_y(:),'ro','MarkerSize',3,'MarkerFace','b');
hold off

