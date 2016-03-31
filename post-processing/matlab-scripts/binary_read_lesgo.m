% Reads in lesgo's binary output data files to facilitate post-processing
%
% author: Joel Bretheim, jbretheim@gmail.com
%
% requires:  lesgo_param.out (in working directory)  
%            binary output files (also, user must specify which snapshot)

clear all; close all; clc;

% specify which files to load
avgVelocities    = 1;
reynoldsStresses = 1;
snapshots        = 1;

% specify file names (must choose a particular velocity snapshot)
str_names_snap = dir('./output/binary_vel.1000.dat.c*');
str_names_avg  = dir('./output/binary_vel_avg.dat*');
str_names_rs   = dir('./output/binary_rs.dat*');

% read in computational domain parameters from lesgo_param.out 
[nx,ny,nz2,nz_tot,nproc,z_i,L_x,L_y,L_z,dx,dy,dz,u_star,nu_molec] = getParams('lesgo_param.out');

% (number of cores used by simulation)
cores = nproc;
dummycores=1:1:cores;
zmin_buf=dummycores*(nz2-1)-nz2+2;
zmax_buf=dummycores*(nz2-1)+1;

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% fetch average velocity fields
if avgVelocities
    [u,v,w] = getAvgVel(nx,ny,nz2,cores,str_names_avg,zmin_buf,zmax_buf);   
end

% fetch instantaneous snapshot velocity fields
if snapshots
    [ubig,vbig,wbig] = getSnap(nx,ny,nz2,cores,str_names_snap,zmin_buf,zmax_buf);   
end

% fetch Reynolds stresses
if reynoldsStresses
    [uu,vv,ww,uw,vw,uv] = getReyStress(nx,ny,nz2,cores,str_names_rs,zmin_buf,zmax_buf);
end

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% build grid
x    = 0.0  : dx : L_x-dx;
y    = 0.0  : dy : L_y-dy;
z_w  = 0.0  : dz : L_z;     % for avg vels, avg (rs's, dudz,dvdz,txz,tyz)
z_uv = dz/2 : dz : L_z;     % for inst vels, avg (txx,txy,tzz,txy)

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% averages across x and y directions (already time-averaged)
if avgVelocities
    uMean = squeeze(mean(mean(u)));
    vMean = squeeze(mean(mean(v)));
    wMean = squeeze(mean(mean(w)));
end

if reynoldsStresses
    uuMean = squeeze(mean(mean(uu)));
    vvMean = squeeze(mean(mean(vv)));
    wwMean = squeeze(mean(mean(ww)));
    uwMean = squeeze(mean(mean(uw)));
    vwMean = squeeze(mean(mean(vw)));
    uvMean = squeeze(mean(mean(uv)));
end

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% basic plots
figure
plot(z_w, uMean,'b')
print('-dpdf','mvp')

figure
kappa = 0.4;  z0 = .0001;
loglaw = 1/kappa * log(z_w ./ z0);    % rough wall
semilogx( z_w, loglaw,'k')
hold on
semilogx(z_w, uMean,'ob')

figure
plot(z_w, uuMean,'ob')

figure
pcolor(ubig(:,:,1))
shading interp; colorbar;









