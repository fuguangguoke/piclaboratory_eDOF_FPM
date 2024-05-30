%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generation of angular spectrum transfer function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H=getH(u,dist, pxsize, wavlen)
[N,M] = size(u);    % size of the wavefront

kx = pi/pxsize*(-1:2/M:1-2/M);
ky = pi/pxsize*(-1:2/N:1-2/N);
[KX,KY] = meshgrid(kx,ky);

k = 2*pi/wavlen;   % wave number
KX_m = KX;
KY_m = KY;
ind = (KX.^2+KY.^2 >= k^2);
KX_m(ind) = 0;
KY_m(ind) = 0;

H = exp(1i*dist*sqrt(k^2-KX_m.^2-KY_m.^2));
end