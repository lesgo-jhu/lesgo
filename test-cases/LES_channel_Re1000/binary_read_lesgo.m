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
domain_snapshots = false;
x_snapshots      = false;
y_snapshots      = false;
z_snapshots      = false;
points           = false;
dns_profiles     = true; % Re_tau = 1000 channel flow

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
if dns_profiles
    dns_data = importdata('dns_profiles.txt');
    dns_data = dns_data.data;
    dns_z  = dns_data(:,1)/1000; % z-plus -> z/h
    dns_u  = dns_data(:,2); % u-plus
    dns_uw = dns_data(:,3);
    dns_uu = dns_data(:,4);
    dns_ww = dns_data(:,5);
    dns_vv = dns_data(:,6);
    dns_tau= dns_data(:,8);
    dns_tot= dns_data(:,9);
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
    
    txzMean = squeeze(mean(mean(txz)));
end

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% basic plots
if avgVelocities
    figure
    %subplot(1,2,1)
    plot(p.z_uv,uMean,'b')
    hold on
    plot(2-p.z_uv,uMean,'r')
    hold off
    ylim([0,30])
    xlabel('z','interpreter','tex')
    ylabel('<u>','interpreter','tex')

    figure
    %subplot(1,2,2)
    kappa = 0.41;  %z0 = .0001186;
    ustar = 0.05; nu = 5e-5; B = 5.2;
    %loglaw = 1/kappa * log(p.z_uv ./ z0);    % rough wall
    loglaw = 1/kappa * log(p.z_uv * ustar / nu) + B; % smooth wall
    semilogx(p.z_uv,loglaw,'k')
    hold on
    semilogx(p.z_uv,uMean,'ob')
    semilogx(dns_z,dns_u,'r')
    %semilogx(2-p.z_uv,uMean,'or')
    hold off
    xlim([0.01,1])
    ylim([0,25])
    xlabel('z','interpreter','tex')
    ylabel('<u>','interpreter','tex')
    legend('Log Law','LES','Location','best')

    figure
    plot(p.z_w, -uwMean,'ob')
    hold on
    %plot(2-p.z_w, uwMean,'or')
    plot(p.z_w, -txzMean,'oc')
    %plot(2-p.z_w, txzMean,'oy')
    plot(p.z_w, -txzMean - uwMean,'ko')
    plot(dns_z,-dns_uw,'b')
    plot(dns_z,dns_tau,'c')
    plot(dns_z,dns_tau-dns_uw,'k')
    plot(p.z_w, (1-p.z_w),'k')
    hold off
    xlabel('z','interpreter','tex')
    ylabel('<u''w''>','interpreter','tex')
    
    figure
    plot(p.z_uv, uuMean, 'ob')
    hold on
    plot(dns_z,dns_uu,'b')
    %plot(2-p.z_uv, uuMean,'or')
    ylim([0,8])
    xlabel('z','interpreter','tex')
    ylabel('<u''u''>','interpreter','tex')
    
    figure
    plot(p.z_uv, vvMean, 'ob')
    hold on
    plot(dns_z,dns_vv, 'b')
    %plot(2-p.z_uv, vvMean,'or')
    ylim([0,4])
    xlabel('z','interpreter','tex')
    ylabel('<v''v''>','interpreter','tex')
    
    figure
    plot(p.z_w, wwMean, 'ob')
    hold on
    plot(dns_z,dns_ww,'b')
    %plot(2-p.z_w, wwMean,'or')
    ylim([0,2])
    xlabel('z','interpreter','tex')
    ylabel('<w''w''>','interpreter','tex')
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