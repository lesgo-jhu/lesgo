function [t,u,v,w] = getPoint(xloc,yloc,zloc)
%GETPOINT reads the time evolution of the velocity at a point
%   [t,u,v,w] = getPoint(xloc,yloc,zloc) reads the time evolution of the velocity 
%   at a point (xloc, yloc, zloc)
    
fname = ['./output/vel.x-',sprintf('%0.5f',xloc),'.y-',sprintf('%0.5f',yloc),...
    '.z-',sprintf('%0.5f',zloc),'.dat'];
dummy=dlmread(fname);
t = dummy(:,1);
u = dummy(:,2);
v = dummy(:,3);
w = dummy(:,4);

end

