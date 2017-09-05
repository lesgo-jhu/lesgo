% Reads in lesgo's binary output data files to facilitate post-processing
%
% author: Joel Bretheim
%         Thanks to Richard Stevens for providing the basic scanning routine
%         which is embedded in the get*.m functions
%
% requires:  lesgo_param.out (in working directory)  
%            binary output files (also, user must specify which snapshot)

clear all; close all; clc;

% specify which files to load
avgVelocities    = true;
dns_profiles     = false; % Re_tau = 1000 channel flow
domain_snapshots = false;
x_snapshots      = false;
y_snapshots      = false;
z_snapshots      = false;
points           = false;

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
    nz = length(uMean);
end

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% basic plots
if avgVelocities
    figure
    plot(p.z_uv(1:nz-1),uMean(1:nz-1),'b')
    hold on
    plot(2-p.z_uv(1:nz-1),-uMean(1:nz-1),'r')
    if dns_profiles
        hold on
        plot(dns_z,dns_u,'r')
    end
    %ylim([-1,1])
    xlabel('z','interpreter','tex')
    ylabel('<u>','interpreter','tex')

    figure
    kappa = 0.41;  nu = 0.001; ustar = sqrt(-txzMean(1)); Re_tau = ustar/nu;
    loglaw = (1/kappa * log(ustar * p.z_uv(1:nz-1) ./ nu) + 5.0);    % rough wall
    semilogx(p.z_uv(1:nz-1)*ustar/nu,loglaw,'k')
    hold on
    semilogx(p.z_uv(1:nz-1)*ustar/nu,p.z_uv(1:nz-1)*ustar/nu,'k--')
    semilogx(p.z_uv(1:nz-1)*ustar/nu,(uMean(1:nz-1)+1)/ustar,'ob')
    semilogx((2-p.z_uv(1:nz-1))*ustar/nu,(1-uMean(1:nz-1))/ustar,'og')
    if dns_profiles
      semilogx(dns_z,dns_u,'r')
    end
    hold off
    xlim([1,2*ustar/nu])
    ylim([0,2/ustar])
    xlabel('z^+','interpreter','tex')
    ylabel('<u>^+','interpreter','tex')
    if dns_profiles
        legend('Log Law','LES','DNS','Location','best')
    else
        legend('Log Law','LES','Location','best')
    end

    figure
    plot(p.z_w*ustar/nu, -uwMean/ustar^2,'ob')
    hold on
    plot(p.z_w*ustar/nu, -txzMean/ustar^2,'oc')
    plot(p.z_w*ustar/nu, (-txzMean - uwMean)/ustar^2,'ko')
    if dns_profiles
        plot(dns_z,-dns_uw,'b')
        plot(dns_z,dns_tau,'c')
        plot(dns_z,dns_tau-dns_uw,'k')
    end
    plot(p.z_w*ustar/nu, ones(nz),'k')
    hold off
    xlabel('z^+','interpreter','tex')
    ylabel('<u''w''>^+','interpreter','tex')
    
    figure
    plot(p.z_uv(1:nz-1)*ustar/nu, uuMean(1:nz-1)/ustar^2, 'ob')
    if dns_profiles
        hold on
        plot(dns_z,dns_uu,'b')
    end
    %ylim([0,8])
    xlabel('z^+','interpreter','tex')
    ylabel('<u''u''>^+','interpreter','tex')
    
    figure
    plot(p.z_uv(1:nz-1)*ustar/nu, vvMean(1:nz-1)/ustar^2, 'ob')
    if dns_profiles
        hold on
        plot(dns_z,dns_vv, 'b')
    end
    %ylim([0,4])
    xlabel('z^+','interpreter','tex')
    ylabel('<v''v''>^+','interpreter','tex')
    
    figure
    plot(p.z_w(1:nz-1)*ustar/nu, wwMean(1:nz-1)/ustar^2, 'ob')
    if dns_profiles
        hold on
        plot(dns_z,dns_ww,'b')
    end
    %ylim([0,2])
    xlabel('z^+','interpreter','tex')
    ylabel('<w''w''>^+','interpreter','tex')
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
