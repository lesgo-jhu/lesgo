function [u2, v2, w2, uw, vw, uv] = getAvgVel2( nx,ny,nz2,nproc,zmin_buf,zmax_buf)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    for i=1:nproc
    
    % Open the file
    fname = ['./output/vel2_avg.c',num2str(i-1),'.bin'];
    fid=fopen(fname,'r');
    if (fid < 0) 
        error('getSnap:fname',['Could not open file ',fname]);
    end

    % Determine the interval of the matrix where the data should be stored
    zmin=zmin_buf(i);
    zmax=zmax_buf(i);
    
    % Scan the data
    dummy=fread(fid,nx*ny*nz2, 'double','s'); 
    u2(1:nx,1:ny,zmin:zmax)=reshape(dummy,nx,ny,nz2);
    dummy=fread(fid,nx*ny*nz2, 'double','s'); 
    v2(1:nx,1:ny,zmin:zmax)=reshape(dummy,nx,ny,nz2);
    dummy=fread(fid,nx*ny*nz2, 'double','s'); 
    w2(1:nx,1:ny,zmin:zmax)=reshape(dummy,nx,ny,nz2);
    dummy=fread(fid,nx*ny*nz2, 'double','s'); 
    uw(1:nx,1:ny,zmin:zmax)=reshape(dummy,nx,ny,nz2);
    dummy=fread(fid,nx*ny*nz2, 'double','s'); 
    vw(1:nx,1:ny,zmin:zmax)=reshape(dummy,nx,ny,nz2);
    dummy=fread(fid,nx*ny*nz2, 'double','s'); 
    uv(1:nx,1:ny,zmin:zmax)=reshape(dummy,nx,ny,nz2);
    
    fclose(fid);
    end
end

