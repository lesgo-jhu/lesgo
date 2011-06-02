clear;
close all;

nproc=8;
zbottom=0.31746;
thresh=1e-9;
sets=[200:200:4000]';
nsets=length(sets);
for ns=1:nsets
  istart=0;
  iend=0;
  for n=1:nproc
    fname=strcat("uvw_avg_z-",int2str(sets(1)),"k-",int2str(sets(ns)),"k.dat.c",int2str(n-1));

    fid=fopen(fname,'r');

    dat = fscanf(fid, '%g %g %g %g ', [4 inf]);
    fclose(fid);
    dat=dat';
    ndat=size(dat,1);
    istart=iend+1;
    iend = istart + ndat - 1 -1; % additional -1 to skip redundant value a top of domain slice
    z(ns,istart:iend) = dat(1:ndat-1,1);
    u(ns,istart:iend) = dat(1:ndat-1,2);
    v(ns,istart:iend) = dat(1:ndat-1,3);
    w(ns,istart:iend) = dat(1:ndat-1,4);
  end
  
end
npoints=size(z,2);
for ns=1:nsets
  for n=1:npoints
    if(u(nsets,n) < thresh)
      u_err(ns,n)=0.;
    else
      u_err(ns,n) = (u(ns,n) - u(nsets,n))/u(nsets,n);
    end

    if(v(nsets,n) < thresh)
      v_err(ns,n)=0.;
    else
      v_err(ns,n) = (v(ns,n) - v(nsets,n))/v(nsets,n);
    end

    if(w(nsets,n) < thresh)
      w_err(ns,n)=0.;
    else
      w_err(ns,n) = (w(ns,n) - w(nsets,n))/w(nsets,n);
    end
  end
end

%  Remove no flow bottom region
z=z-zbottom;

for ns=1:nsets
  fname=strcat("uvw_avg_z-",int2str(sets(1)),"k-",int2str(sets(ns)),"k.dat");
  fid=fopen(fname,'w');
  for n=1:npoints
    if(z(ns,n) >= 0)
      fprintf(fid,'%g \t %g \t %g \t %g \t %g \t %g \t %g\n', ...
        z(ns,n), u(ns,n), v(ns,n), w(ns,n), u_err(ns,n), v_err(ns,n), w_err(ns,n));
    end
  end
  fclose(fid);  
end

%  Compute Reynolds stress data

for ns=1:nsets
  istart=0;
  iend=0;
  for n=1:nproc
    fname=strcat("rs_z-",int2str(sets(1)),"k-",int2str(sets(ns)),"k.dat.c",int2str(n-1));

    fid=fopen(fname,'r');

    dat = fscanf(fid, '%g %g %g %g %g %g %g', [7 inf]);
    fclose(fid);
    dat=dat';
    ndat=size(dat,1);
    istart=iend+1;
    iend = istart + ndat - 1 -1; % additional -1 to skip redundant value a top of domain slice
    z(ns,istart:iend) = dat(1:ndat-1,1);
    up2(ns,istart:iend) = dat(1:ndat-1,2);
    vp2(ns,istart:iend) = dat(1:ndat-1,3);
    wp2(ns,istart:iend) = dat(1:ndat-1,4);
    upwp(ns,istart:iend) = dat(1:ndat-1,5);
    vpwp(ns,istart:iend) = dat(1:ndat-1,6);
    upvp(ns,istart:iend) = dat(1:ndat-1,7);

  end
  
end

npoints=size(z,2);
for ns=1:nsets
  for n=1:npoints
    if(up2(nsets,n) < thresh)
      up2_err(ns,n)=0.;
    else
      up2_err(ns,n) = (up2(ns,n) - up2(nsets,n))/up2(nsets,n);
    end

    if(vp2(nsets,n) < thresh)
      vp2_err(ns,n)=0.;
    else
      vp2_err(ns,n) = (vp2(ns,n) - vp2(nsets,n))/vp2(nsets,n);
    end

    if(wp2(nsets,n) < thresh)
      wp2_err(ns,n)=0.;
    else
      wp2_err(ns,n) = (wp2(ns,n) - wp2(nsets,n))/wp2(nsets,n);
    end

    if(upwp(nsets,n) < thresh)
      upwp_err(ns,n)=0.;
    else
      upwp_err(ns,n) = (upwp(ns,n) - upwp(nsets,n))/upwp(nsets,n);
    end

    if(vpwp(nsets,n) < thresh)
      vpwp_err(ns,n)=0.;
    else
      vpwp_err(ns,n) = (vpwp(ns,n) - vpwp(nsets,n))/vpwp(nsets,n);
    end

    if(upvp(nsets,n) < thresh)
      upvp_err(ns,n)=0.;
    else
      upvp_err(ns,n) = (upvp(ns,n) - upvp(nsets,n))/upvp(nsets,n);
    end

  end
end

%  Remove no flow bottom region
z=z-zbottom;

for ns=1:nsets
  fname=strcat("rs_z-",int2str(sets(1)),"k-",int2str(sets(ns)),"k.dat");
  fid=fopen(fname,'w');
  for n=1:npoints
    if(z(ns,n) >= 0)
      fprintf(fid,'%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n', ...
        z(ns,n), up2(ns,n), vp2(ns,n), wp2(ns,n), upwp(ns,n), vpwp(ns,n), upvp(ns,n), ...
        up2_err(ns,n), vp2_err(ns,n), wp2_err(ns,n), upwp_err(ns,n), vpwp_err(ns,n), upvp_err(ns,n));
    end
  end
  fclose(fid);  
end

