%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% complex Hessian denoiser for pupil function enhancement
% input: 
% G: noisy complex pupil function
% tau: threshold value for denoising
% num: number of Hessian denoiser iteration
% output: 
% x: denoised complex pupil function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x]=HessianOpt(G,tau,num)
[M,N]=size(G);
hxx=(psf2otf([1 -2 1],[M N]));
hxy=(psf2otf([-1 1;1 -1],[M N]));
hyy=(psf2otf([1; -2; 1],[M N]));
hyx=(psf2otf([1 -1; -1 1],[M N]));
s=fft2(G);
W_xx=ifft2(s.*hxx);
W_xy=ifft2(s.*hxy);
W_yy=ifft2(s.*hyy);
W_yx=ifft2(s.*hyx);


y=zeros([M N]);
v_xx=zeros([M N]);
v_xy=zeros([M N]);
v_yy=zeros([M N]);
v_yx=zeros([M N]);



mu=1;eta=1;alpha=1;beta=1;
for k=1:num
    x=1/(1+mu)*(G+mu*ifft2(s+1/mu*y));
    numer1=(fft2(x)-1/mu*y);
    numer2=conj(hxx).*fft2(W_xx+1/eta*v_xx)+conj(hxy).*fft2(W_xy+1/eta*v_xy)+conj(hyy).*fft2(W_yy+1/eta*v_yy)+conj(hyx).*fft2(W_yx+1/eta*v_yx);
    numer_s=mu*numer1+eta*numer2;
    demoni_s=mu+eta*(hxx.*conj(hxx)+hxy.*conj(hxy)+hyy.*conj(hyy)+hyx.*conj(hyx));
    s=numer_s./demoni_s;
    
    Pxx=ifft2(hxx.*s)-1/eta*v_xx;Pxy=ifft2(hxy.*s)-1/eta*v_xy;
    Pyy=ifft2(hyy.*s)-1/eta*v_yy;Pyx=ifft2(hyx.*s)-1/eta*v_yx;
    [W_xx,W_xy,W_yy,W_yx]=softshrink(Pxx,Pxy,Pyy,Pyx,4*tau/eta,'isotropic');
    
    y=y+alpha*mu*(s-fft2(x));
    v_xx=v_xx+beta*eta*(W_xx-ifft2(hxx.*s));
    v_xy=v_xy+beta*eta*(W_xy-ifft2(hxy.*s));
    v_yy=v_yy+beta*eta*(W_yy-ifft2(hyy.*s));
    v_yx=v_yx+beta*eta*(W_xx-ifft2(hyx.*s));
end
end
function [W1,W2,W3,W4]=softshrink(x1,x2,x3,x4,beta,methodTV)
switch methodTV
    case 'isotropic'
        s = sqrt(abs(x1).^2 + abs(x2).^2+abs(x3).^2+abs(x4).^2);
        W1 = (x1./(s + eps)).*max(s - beta,0);
        W2 = (x2./(s + eps)).*max(s - beta,0);
        W3 = (x3./(s + eps)).*max(s - beta,0);
        W4 = (x4./(s + eps)).*max(s - beta,0);
    case 'anisotropic'
        W1 = sign(x1).*max(abs(x1) - beta,0);
        W2 = sign(x2).*max(abs(x2) - beta,0);
        W3 = sign(x3).*max(abs(x3) - beta,0);
        W4 = sign(x4).*max(abs(x4) - beta,0);
end
end
