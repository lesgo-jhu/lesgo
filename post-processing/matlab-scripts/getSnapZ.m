function [ u,v,w ] = getSnapZ( step,loc,nx,ny,nproc)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    for i=1:nproc
    
    % Open the file
    fname = ['./output/vel.z-',sprintf('%0.5f',loc),'.',num2str(step),'.bin'];
    fid=fopen(fname,'r');
    if (fid < 0) 
        error('getSnap:fname',['Could not open file ',fname]);
    end

    % Scan the data
    dummy=fread(fid,nx*ny, 'double','s');
    u=reshape(dummy,nx,ny);
    dummy=fread(fid,nx*ny, 'double','s'); 
    v=reshape(dummy,nx,ny);
    dummy=fread(fid,nx*ny, 'double','s'); 
    w=reshape(dummy,nx,ny);
    
    fclose(fid);

    end

end

