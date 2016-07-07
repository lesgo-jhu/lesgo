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
p = getParams('output/lesgo_param.out');

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% fetch average velocity fields
if avgVelocities
    [u,v,w_uv] = getAvgVelUV(p);
    w = getAvgVelW(p);
    [uu,vv,ww,uw,vw,uv] = getReyStress(p);
    [u2,v2,w2,uw2,vw2,uv2] = getAvgVel2(p);
    [txx,txy,tyy,txz,tyz,tzz] = getAvgTau(p);
    [fx,fy,fz] = getAvgForce(p);
    CS = getCSOpt2(p);
end

% fetch instantaneous snapshot velocity fields
if domain_snapshots
    [ubig,vbig,wbig] = getSnap(p,snap_time);
end
if x_snapshots
    [ux,vx,wx] = getSnapX(p,snap_time,xloc);
end
if y_snapshots
    [uy,vy,wy] = getSnapY(p,snap_time,yloc);
end
if z_snapshots
    [uz,vz,wz] = getSnapZ(p,snap_time,zloc);
end
if points
    [t,up,vp,wp] = getPoint(ploc(1),ploc(2),ploc(3));
end

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% averages across p.x and p.y directions (already time-averaged)
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
    plot(p.z_w,uMean,'b')
    xlabel('z','interpreter','tex')
    ylabel('<u>','interpreter','tex')

    subplot(1,2,2)
    kappa = 0.4;  z0 = .0001;
    loglaw = 1/kappa * log(p.z_w ./ z0);    % rough wall
    semilogx(p.z_w,loglaw,'k')
    hold on
    semilogx(p.z_w,uMean,'ob')
    hold off
    xlabel('z','interpreter','tex')
    ylabel('<u>','interpreter','tex')
    legend('Log Law','LES')

    figure
    plot(p.z_uv, uuMean,'ob')
    xlabel('z','interpreter','tex')
    ylabel('<u''u''>','interpreter','tex')
end

if domain_snapshots
    figure
    [X,Y] = meshgrid(p.x,p.y);
    pcolor(X,Y,ubig(:,:,4)')
    xlabel('x')
    ylabel('y')
    title(['Streamwise velocity at z = ',num2str(p.z_uv(4))])
    shading interp; colorbar;
end
    
if x_snapshots
    figure
    [Y,Z] = meshgrid(p.y,p.z_uv);
    pcolor(Y,Z,ux')
    xlabel('y')
    ylabel('z')
    title(['Streamwise velocity at x = ',num2str(xloc)])
    shading interp; colorbar;
end

if y_snapshots
    figure
    [X,Z] = meshgrid(p.x,p.z_uv);
    pcolor(X,Z,uy')
    xlabel('x')
    ylabel('z')
    title(['Streamwise velocity at y = ',num2str(yloc)])
    shading interp; colorbar;
end

if z_snapshots
    figure
    [X,Y] = meshgrid(p.x,p.y);
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