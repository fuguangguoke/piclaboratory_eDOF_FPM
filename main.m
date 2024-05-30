clear;close all;clc
addpath('./function/');
addpath('./data/');
%% parameters
load('aberration.mat')
Pupil=Pupil0;%                       system pupil function
lamuda=0.632;%                       wavelength,unit is micron
k_lamuda=2*pi/lamuda;%               wave vector
D_led=4*1000;%                       interval between LEDs
h=86*1000;%                          height from LED to sample
ledMM=9;ledNN=9;%                  scanning array
ledM=ledMM;ledN=ledNN;

pixel_size=6.5;%                     pixel size of camera
mag=4;%                              magnification of objective
NA=0.1;%                             numerical aperture of objective
M=128;N=128;%                        segment size
D_pixel=pixel_size/mag;%             pixel size of sample
kmax=NA*k_lamuda;


dist=200;%                           200 um defocus
Niter1=100;%                         iteration times of FPM

Rcam=lamuda/NA*mag/2/pixel_size;
RLED=NA*sqrt(D_led^2+h^2)/D_led;
Roverlap=1/pi*(2*acos(1/2/RLED)-1/RLED*sqrt(1-(1/2/RLED)^2));

maxr_ring=7;
MAGimg=ceil(1+2*D_pixel*maxr_ring*D_led/sqrt((maxr_ring*D_led)^2+h^2)/lamuda);% factor of upsampling
MM=M*MAGimg;NN=N*MAGimg;


% segment position
x=0;
objdx=x*D_pixel;
y=0;
objdy=y*D_pixel;

%% generate object
I0=im2double(imread('I0.bmp'));

P0=im2double((imread('P0.bmp')));
I0=imresize(I0,[MM,NN]);
P0=imresize(P0,[MM,NN]);
o1=sqrt(I0).*exp(sqrt(-1).*P0);
O1=fftshift(fft2(o1));

%% generate coordinate of frequency domain
[Fx1,Fy1]=meshgrid(-(N/2):(N/2-1),-(M/2):(M/2-1));
Fx1=Fx1./(N*D_pixel).*(2*pi);
Fy1=Fy1./(M*D_pixel).*(2*pi);
Fx2=Fx1.*Fx1;
Fy2=Fy1.*Fy1;
Fxy2=Fx2+Fy2;
Pupil00=zeros(M,N);
Pupil00(Fxy2<=(kmax^2))=1;
[Fxx1,Fyy1]=meshgrid(-(NN/2):(NN/2-1),-(MM/2):(MM/2-1));
Fxx1=Fxx1(1,:)./(N*D_pixel).*(2*pi);
Fyy1=Fyy1(:,1)./(M*D_pixel).*(2*pi);

%% LED position in K space
lit_cenv = 4;
lit_cenh = 4;
vled = (0:8)-lit_cenv;
hled = (0:8)-lit_cenh;
[hhled,vvled] = meshgrid(hled,vled);

k=zeros(ledMM,ledNN);
for i=1:ledMM
    for j=1:ledNN
        k(i,j)=j+(i-1)*ledNN;
    end
end

v = (-hhled*D_led-objdx)./sqrt((-hhled*D_led-objdx).^2+(vvled*D_led-objdy).^2+h.^2);
u = (vvled*D_led-objdy)./sqrt((-hhled*D_led-objdx).^2+(vvled*D_led-objdy).^2+h.^2);
NAillu=sqrt(u.^2+v.^2);  
ledpos_true=zeros(ledMM,ledNN,2);

for i=1:ledMM
    for j=1:ledNN
        if k(i,j)~=0
        Fx1_temp=abs(Fxx1-k_lamuda*u(i,j));
        ledpos_true(i,j,1)=find(Fx1_temp==min(Fx1_temp));
        Fy1_temp=abs(Fyy1-k_lamuda*v(i,j));
        ledpos_true(i,j,2)=find(Fy1_temp==min(Fy1_temp));
        end
    end
end

%% image recording
Isum=zeros(M,N,ledMM*ledNN);
H=getH(zeros([M N]),dist,D_pixel,lamuda);
for i=1:ledMM
    for j=1:ledNN
        if k(i,j)~=0
            uo=ledpos_true(i,j,1);
            vo=ledpos_true(i,j,2);
            O1P0=O1((vo-M/2):(vo-1+M/2),(uo-N/2):(uo-1+N/2))./(MAGimg^2).*Pupil.*H;
            o10=ifft2(fftshift(O1P0));
            oI0=abs(o10).^2;
            Isum(:,:,k(i,j))=oI0;
        end
    end
end
I_mid=Isum(:,:,(ledMM*ledNN+1)/2);%center image with normal incidence


%% update order
Nled=ledM*ledN;
ord_ijsum=zeros(ledMM,ledNN);
ord_isum=zeros(1,Nled);
ord_jsum=zeros(1,Nled);
ord_ii=(ledMM+1)/2; 
ord_jj=(ledNN+1)/2; 
ord_isum(1,1)=ord_ii;
ord_jsum(1,1)=ord_jj;
ord_ijsum(ord_ii,ord_jj)=1;
led_num=1;
direction=0;

while (min(min(ord_ijsum))==0)
    led_num=led_num+1;
    direction2=direction+1;
    ord_ii2=round(ord_ii+sin(pi/2*direction2));
    ord_jj2=round(ord_jj+cos(pi/2*direction2));
    if (ord_ijsum(ord_ii2,ord_jj2)==0)
        direction=direction2;
    end
    ord_ii=round(ord_ii+sin(pi/2*direction));
    ord_jj=round(ord_jj+cos(pi/2*direction));
    ord_isum(1,led_num)=ord_ii;
    ord_jsum(1,led_num)=ord_jj;
    ord_ijsum(ord_ii,ord_jj)=1;
end

A_GT=sqrt(I0);
FF=fft2(A_GT);
[Row,Col]=size(FF);
%% initial setup

S=ones(MM,NN);
s=fftshift(fft2(S));
Pupil=Pupil00;

% initial parameters
alpha=1;
beta=1;
gamma=0.3;
eta=1;
taup=5e-4;

% initial value for w-update
w=cell(Nled,1);
for led_num=1:1:Nled
    w{led_num}=complex(zeros(M,N)+j*zeros(M,N));
end

% initial value for q-update
q=cell(Nled,1);
for led_num=1:1:Nled
    q{led_num}=complex(zeros(M,N)+j*zeros(M,N));
end

r=0; 
s_store=cell(Niter1,1);
w_store=cell(Niter1,1);


for iter=1:Niter1
    error_now=0;
    q1_sum=complex(zeros(MM,NN)+j*zeros(MM,NN));
    q2_sum=complex(zeros(M,N)+j*zeros(M,N));
    p_sum=complex(zeros(MM,NN)+j*zeros(MM,NN));
    s_sum=complex(zeros(M,N)+j*zeros(M,N));
   %% ADMM-FPM update
   % q-update
   for led_num=1:1:Nled 
                i=ord_isum(led_num); j=ord_jsum(led_num);
                uo=ledpos_true(i,j,1); vo=ledpos_true(i,j,2);
                I_sqrt=sqrt(Isum(:,:,k(i,j)));  
                q{led_num}=s((vo-M/2):(vo-1+M/2),(uo-N/2):(uo-1+N/2)).*Pupil./(MAGimg^2)-w{led_num};
                Q=ifft2(fftshift(q{led_num}));              
                q{led_num}= fftshift(fft2((I_sqrt+alpha*abs(Q))/(1+alpha).*sign(Q)));
                if iter>2 
                r=r+(norm(reshape(q{led_num}-s_store{iter-1}((vo-M/2):(vo-1+M/2),(uo-N/2):(uo-1+N/2)).*Pupil./(MAGimg^2)+w_store{iter-1}{led_num},1,[]),2)).^2;
                end
    end
                          
    % s-update
    for led_num=1:1:Nled   
                i=ord_isum(led_num); j=ord_jsum(led_num);
                uo=ledpos_true(i,j,1); vo=ledpos_true(i,j,2);
                q1_sum((vo-M/2):(vo-1+M/2),(uo-N/2):(uo-1+N/2))=q1_sum((vo-M/2):(vo-1+M/2),(uo-N/2):(uo-1+N/2))+conj(Pupil).*(q{led_num}+w{led_num}).*(MAGimg^2);
                p_sum((vo-M/2):(vo-1+M/2),(uo-N/2):(uo-1+N/2))=p_sum((vo-M/2):(vo-1+M/2),(uo-N/2):(uo-1+N/2))+abs(Pupil).^2;
    end
    s=q1_sum./(gamma/alpha + p_sum);
    s_store{iter}=s;
    
    % p-update
    for led_num=1:1:Nled 
                i=ord_isum(led_num); j=ord_jsum(led_num);
                uo=ledpos_true(i,j,1); vo=ledpos_true(i,j,2);
                q2_sum=q2_sum + conj(s((vo-M/2):(vo-1+M/2),(uo-N/2):(uo-1+N/2)) ./(MAGimg^2)).*(q{led_num}+w{led_num});
                s_sum=s_sum+ abs( s((vo-M/2):(vo-1+M/2),(uo-N/2):(uo-1+N/2))).^2;

    end
    Pupil=q2_sum./(beta+ s_sum );
   
    %% Hessian denoiser for pupil
    Pupil=HessianOpt(Pupil,taup,5);

    
    % support constraint for pupil
    Pupil=Pupil00.*exp(1i*angle(Pupil));
    
    % w-update
    for led_num=1:1:Nled 
                i=ord_isum(led_num); j=ord_jsum(led_num);
                uo=ledpos_true(i,j,1); vo=ledpos_true(i,j,2);
                w{led_num} = w{led_num} + eta * (q{led_num}-s((vo-M/2):(vo-1+M/2),(uo-N/2):(uo-1+N/2)).*Pupil./(MAGimg^2));     
    end
    w_store{iter}=w;
    
    disp([num2str(iter),'/',num2str(Niter1)]);
    
  
    

end
  % reconstruction results
  figure;
  nexttile;imshow(abs(o1),[]);title('Object Amplitude ground truth')
  nexttile;imshow(angle(o1),[]);title('Object Phase ground truth')
  nexttile;imshow(I_mid,[]);title([' BF image of ',num2str(dist),' um defocus'])
  s_A=abs(ifft2(fftshift(s)));
  s_p=angle(ifft2(fftshift(s)))+0.5;
  nexttile;imshow(s_A,[]);
  title(['Iteration No. = ',int2str(iter),' recovered Object Amplitude']);
  nexttile;imshow(s_p,[]);colorbar
  title(['Iteration No. = ',int2str(iter),' recovered Object Phase']);
  nexttile;imshow(angle(Pupil),[-pi pi]);
  title(['Iteration No. = ',int2str(iter),' recovered Pupil Phase']);
