% Reads in lesgo's binary output data files to facilitate post-processing
%
% author: Joel Bretheim, jbretheim@gmail.com
%         Thanks to Richard Stevens for providing the basic scanning routine
%         which is embedded in the get*.m functions
%
% requires:  lesgo_param.out (in working directory)  
%            binary output files (also, user must specify which snapshot)

clear all; close all; clc;

% specify which files to load
avgVelocities    = true;
domain_snapshots = true;
x_snapshots      = true;
y_snapshots      = true;
z_snapshots      = true;
points           = true;

% specify file names (must choose a particular velocity snapshot)
snap_time = 1;
xloc = 0.1;
yloc = 0.1;
zloc = 0.1;
ploc = [0.1,0.1,0.1];

% read in computational domain parameters from lesgo_param.out 
[nx,ny,nz2,nz_tot,nproc,z_i,L_x,L_y,L_z,dx,dy,dz,u_star,nu_molec] = getParams('lesgo_param.out');

% (number of nproc used by simulation)
dummynproc=1:1:nproc;
zmin_buf=dummynproc*(nz2-1)-nz2+2;
zmax_buf=dummynproc*(nz2-1)+1;

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% fetch average velocity fields
if avgVelocities
    [u,v,w_uv] = getAvgVelUV(nx,ny,nz2,nproc,zmin_buf,zmax_buf);
    w = getAvgVelW(nx,ny,nz2,nproc,zmin_buf,zmax_buf);
    [uu,vv,ww,uw,vw,uv] = getReyStress(nx,ny,nz2,nproc,zmin_buf,zmax_buf);
    [u2,v2,w2,uw2,vw2,uv2] = getAvgVel2(nx,ny,nz2,nproc,zmin_buf,zmax_buf);
    [txx,txy,tyy,txz,tyz,tzz] = getAvgTau(nx,ny,nz2,nproc,zmin_buf,zmax_buf);
    fx = getAvgForce(nx,ny,nz2,nproc,zmin_buf,zmax_buf);
    CS = getCSOpt2(nx,ny,nz2,nproc,zmin_buf,zmax_buf);
end

% fetch instantaneous snapshot velocity fields
if domain_snapshots
    [ubig,vbig,wbig] = getSnap(snap_time,nx,ny,nz2,nproc,zmin_buf,zmax_buf);   
end
if x_snapshots
    [ux,vx,wx] = getSnapX(snap_time,xloc,ny,nz2,nproc,zmin_buf,zmax_buf);   
end
if y_snapshots
    [uy,vy,wy] = getSnapY(snap_time,yloc,ny,nz2,nproc,zmin_buf,zmax_buf);   
end
if z_snapshots
    [uz,vz,wz] = getSnapZ(snap_time,zloc,nx,ny,nproc);
end
if points
    [t,up,vp,wp] = getPoint(ploc(1),ploc(2),ploc(3));
end

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% build grid
x    = 0.0  : dx : L_x-dx;
y    = 0.0  : dy : L_y-dy;
z_w  = 0.0  : dz : L_z;     % for avg vels, avg (rs's, dudz,dvdz,txz,tyz)
z_uv = dz/2 : dz : L_z+dz/2;     % for inst vels, avg (txx,txy,tzz,txy)

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% averages across x and y directions (already time-averaged)
if avgVelocities
    uMean = squeeze(mean(mean(u)));
    vMean = squeeze(mean(mean(v)));
    wMean = squeeze(mean(mean(w)));

    uuMean = squeeze(mean(mean(uu)));
    vvMean = squeeze(mean(mean(vv)));
    wwMean = squeeze(mean(mean(ww)));
    uwMean = squeeze(mean(mean(uw)));
    vwMean = squeeze(mean(mean(vw)));
    uvMean = squeeze(mean(mean(uv)));
end

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% basic plots
if avgVelocities
    figure
    subplot(1,2,1)
    plot(z_w,uMean,'b')
    xlabel('z','interpreter','tex')
    ylabel('<u>','interpreter','tex')

    subplot(1,2,2)
    kappa = 0.4;  z0 = .0001;
    loglaw = 1/kappa * log(z_w ./ z0);    % rough wall
    semilogx(z_w,loglaw,'k')
    hold on
    semilogx(z_w,uMean,'ob')
    hold off
    xlabel('z','interpreter','tex')
    ylabel('<u>','interpreter','tex')
    legend('Log Law','LES')

    figure
    plot(z_uv, uuMean,'ob')
    xlabel('z','interpreter','tex')
    ylabel('<u''u''>','interpreter','tex')
end

if domain_snapshots
    figure
    [X,Y] = meshgrid(x,y);
    pcolor(X,Y,ubig(:,:,4)')
    xlabel('x')
    ylabel('y')
    title(['Streamwise velocity at z = ',num2str(z_uv(4))])
    shading interp; colorbar;
end
    
if x_snapshots
    figure
    [Y,Z] = meshgrid(y,z_uv);
    pcolor(Y,Z,ux')
    xlabel('y')
    ylabel('z')
    title(['Streamwise velocity at x = ',num2str(xloc)])
    shading interp; colorbar;
end

if y_snapshots
    figure
    [X,Z] = meshgrid(x,z_uv);
    pcolor(X,Z,uy')
    xlabel('x')
    ylabel('z')
    title(['Streamwise velocity at y = ',num2str(yloc)])
    shading interp; colorbar;
end

if z_snapshots
    figure
    [X,Y] = meshgrid(x,y);
    pcolor(X,Y,uz')
    xlabel('x')
    ylabel('y')
    title(['Streamwise velocity at z = ',num2str(zloc)])
    shading interp; colorbar;
end

if points
    figure
    plot(t,[up,vp,wp])
    xlabel('t')
    ylabel('Velocity')
    legend('u','v','w')
end







