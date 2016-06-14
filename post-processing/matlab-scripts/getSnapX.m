function [ u,v,w ] = getSnapX( step,loc,ny,nz2,nproc,zmin_buf,zmax_buf)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    for i=1:nproc
    
    % Open the file
    fname = ['./output/vel.x-',sprintf('%0.5f',loc),'.',num2str(step),'.c',num2str(i-1),'.bin'];
    fid=fopen(fname,'r');
    if (fid < 0) 
        error('getSnap:fname',['Could not open file ',fname]);
    end

    % Determine the interval of the matrix where the data should be stored
    zmin=zmin_buf(i);
    zmax=zmax_buf(i);

    % Scan the data
    dummy=fread(fid,ny*nz2, 'double','s');
    u(1:ny,zmin:zmax)=reshape(dummy,ny,nz2);
    dummy=fread(fid,ny*nz2, 'double','s'); 
    v(1:ny,zmin:zmax)=reshape(dummy,ny,nz2);
    dummy=fread(fid,ny*nz2, 'double','s'); 
    w(1:ny,zmin:zmax)=reshape(dummy,ny,nz2);
    
    fclose(fid);

    end

end

