function [ u,v,w ] = getSnapX(p,step,loc)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for i=1:p.nproc
    
    % Open the file
    fname = ['./output/vel.x-',sprintf('%0.5f',loc),'.',num2str(step),'.c',num2str(i-1),'.bin'];
    fid=fopen(fname,'r');
    if (fid < 0) 
        error('getSnap:fname',['Could not open file ',fname]);
    end

    % Determine the interval of the matrix where the data should be stored
    zmin=p.zmin_buf(i);
    zmax=p.zmax_buf(i);

    % Scan the data
    dummy=fread(fid,p.ny*p.nz2,'double',p.fmt);
    u(1:p.ny,zmin:zmax)=reshape(dummy,p.ny,p.nz2);
    dummy=fread(fid,p.ny*p.nz2,'double',p.fmt); 
    v(1:p.ny,zmin:zmax)=reshape(dummy,p.ny,p.nz2);
    dummy=fread(fid,p.ny*p.nz2,'double',p.fmt); 
    w(1:p.ny,zmin:zmax)=reshape(dummy,p.ny,p.nz2);
    
    fclose(fid);

end

end

