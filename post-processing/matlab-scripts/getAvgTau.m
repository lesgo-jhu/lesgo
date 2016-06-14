function [txx,txy,tyy,txz,tyz,tzz] = getAvgTau( nx,ny,nz2,nproc,zmin_buf,zmax_buf)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    for i=1:nproc
    
    % Open the file
    fname = ['./output/tau_avg.c',num2str(i-1),'.bin'];
    fid=fopen(fname,'r');
    if (fid < 0) 
        error('getSnap:fname',['Could not open file ',fname]);
    end

    % Determine the interval of the matrix where the data should be stored
    zmin=zmin_buf(i);
    zmax=zmax_buf(i);
    
    % Scan the data
    dummy=fread(fid,nx*ny*nz2, 'double','s'); 
    txx(1:nx,1:ny,zmin:zmax)=reshape(dummy,nx,ny,nz2);
    dummy=fread(fid,nx*ny*nz2, 'double','s'); 
    txy(1:nx,1:ny,zmin:zmax)=reshape(dummy,nx,ny,nz2);
    dummy=fread(fid,nx*ny*nz2, 'double','s'); 
    tyy(1:nx,1:ny,zmin:zmax)=reshape(dummy,nx,ny,nz2);
    dummy=fread(fid,nx*ny*nz2, 'double','s'); 
    txz(1:nx,1:ny,zmin:zmax)=reshape(dummy,nx,ny,nz2);
    dummy=fread(fid,nx*ny*nz2, 'double','s'); 
    tyz(1:nx,1:ny,zmin:zmax)=reshape(dummy,nx,ny,nz2);
    dummy=fread(fid,nx*ny*nz2, 'double','s'); 
    tzz(1:nx,1:ny,zmin:zmax)=reshape(dummy,nx,ny,nz2);
    
    fclose(fid);
    end

end

