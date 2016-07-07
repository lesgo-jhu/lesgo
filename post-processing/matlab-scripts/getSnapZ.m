function [ u,v,w ] = getSnapZ(p,step,loc)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for i=1:p.nproc
    
    % Open the file
    fname = ['./output/vel.z-',sprintf('%0.5f',loc),'.',num2str(step),'.bin'];
    fid=fopen(fname,'r');
    if (fid < 0) 
        error('getSnap:fname',['Could not open file ',fname]);
    end

    % Scan the data
    dummy=fread(fid,p.nx*p.ny, 'double',p.fmt);
    u=reshape(dummy,p.nx,p.ny);
    dummy=fread(fid,p.nx*p.ny, 'double',p.fmt); 
    v=reshape(dummy,p.nx,p.ny);
    dummy=fread(fid,p.nx*p.ny, 'double',p.fmt); 
    w=reshape(dummy,p.nx,p.ny);
    
    fclose(fid);

end

end

