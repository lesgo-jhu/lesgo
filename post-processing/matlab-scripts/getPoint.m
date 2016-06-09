function [t,u,v,w] = getPoint(xloc,yloc,zloc)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    
    fname = ['./output/vel.x-',sprintf('%0.5f',xloc),'.y-',sprintf('%0.5f',yloc),...
        '.z-',sprintf('%0.5f',zloc),'.dat'];
    dummy=dlmread(fname);
    t = dummy(:,1);
    u = dummy(:,2);
    v = dummy(:,3);
    w = dummy(:,4);

end

