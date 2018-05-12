function [u,v,w,vortx,vorty,vortz,p] = getSnap(p,step)
%GETSNAP reads whole domain snapshot
%   [u,v,w,vortx,vorty,vortz,p] = getSnap(p,step) reads a domain snapshot at timestep step. 
%   lesgo parameters are provided as struct p read using p = getParams(...)

% size of record
N = p.nx*p.ny*p.nz_tot;

% Velocity
fname = ['./output/vel.',num2str(step),'.bin'];
fid = fopen(fname,'r');
if (fid < 0) 
    error(['getSnap: Could not open file ',fname]);
end
u = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
v = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
w = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
fclose(fid);

% Vorticity
fname = ['./output/vort.',num2str(step),'.bin'];
fid = fopen(fname,'r');
if (fid < 0) 
    error(['getSnap: Could not open file ',fname]);
end
vortx = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
vorty = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
vortz = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
fclose(fid);

% Pressure
fname = ['./output/pres.',num2str(step),'.bin'];
fid = fopen(fname,'r');
if (fid < 0) 
    error(['getSnap: Could not open file ',fname]);
end
p = reshape(fread(fid,N,'double',p.fmt),p.nx,p.ny,p.nz_tot);
fclose(fid);

end

