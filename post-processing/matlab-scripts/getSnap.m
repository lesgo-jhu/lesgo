function [ ubig,vbig,wbig ] = getSnap( nx,ny,nz2,cores,str,zmin_buf,zmax_buf)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    for i=1:cores

    % Access the different files    
    name2 = str(i).name;
    fid=fopen(strcat('./output/',name2),'r');

    % Determine the core number of the file in order to store the data in the right location
    name2 = regexp(name2, '.c', 'split');  % split of the last numbers behind the core number. This indicates the domain 
    domain=char(name2(2));                % convert the cell structure to a string structure
    domain=str2double(domain);            % convert the string to a numerical value that can be used later on

    % Determine the interval of the matrix where the data should be stored
    zmin=zmin_buf(domain+1);
    zmax=zmax_buf(domain+1);

    % Scan the data
    dummy=fread(fid,nx*ny*nz2, 'double','s');
    ubig(1:nx,1:ny,zmin:zmax)=reshape(dummy,nx,ny,nz2);
    dummy=fread(fid,nx*ny*nz2, 'double','s'); 
    vbig(1:nx,1:ny,zmin:zmax)=reshape(dummy,nx,ny,nz2);
    dummy=fread(fid,nx*ny*nz2, 'double','s'); 
    wbig(1:nx,1:ny,zmin:zmax)=reshape(dummy,nx,ny,nz2);
    
    fclose(fid);

    end
end

