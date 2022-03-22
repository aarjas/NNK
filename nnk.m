function [latcoords,axcoords,elevcoords] = nnk(images,RFdata,nn,showplot)
net3=nn;
data = RFdata;
ns = 15;
nx = 1001;
nz = 601;
DynRan = 30;
deltat = 100/size(images,3);
states = 8;
q = 0.005;
q1 = q^2;
q2 = q^2;
q3 = q^2; 
q4 = q^2; 
H = [[ones(ns,1);zeros(ns+2,1)] [zeros(ns,1);ones(ns,1);0;0] [zeros(2*ns,1);1;0] [zeros(2*ns+1,1);1] zeros(2*ns+2,4)];
G1 = 0.5*deltat^2*eye(4);
G2 = eye(4);
G = [G1;G2];
Q = G*diag([q1 q2 q3 q4])*G';
A = eye(states) + diag(ones(4,1),4);
mk = zeros(8,1);
Pk = 15*eye(8);
IMAGE_env = abs(hilbert(images(:,:,20)'))';
IMAGE_log = single(20*log10(IMAGE_env'));
mc = max(max(IMAGE_log));
IMAGE_env = abs(hilbert(data(:,:,20)'))';
IMAGE_log = single(20*log10(IMAGE_env'));
mc2 = max(max(IMAGE_log));
ym=linspace(0,19.1975,4800);
xm=1:64;
zpos=zeros(1,64);
xpos=linspace(-0.0125,0.0125,64);
hydrophone=1e-3*[0 0.4 0];
zimg=linspace(0,15,601)*1e-3;
ximg=linspace(-12.5,12.5,1001)*1e-3;


  

latcoords = zeros(size(images,3),1);
axcoords = zeros(size(images,3),1);
elevcoords = zeros(size(images,3),1);

if showplot
    figure('Renderer', 'painters', 'Position', [100 100 1001 300],'Color','White')
end

pos=[0 7.5];
for k = 1:size(images,3)
    y = images(:,:,k);
    d = data(:,:,k);
    standd=(d(:)-mean(d(:)))/std(d(:));
    I=mean(abs(d(abs(standd)>norminv(0.99))));
    I=log(I);
    
    [~,idx] = sort(abs(y(:)),'descend');
    [row,col] = ind2sub(size(y),idx(1:ns));
    xx = -12.5+(25/(nx-1))*(row+1);
    zz = (15/(nz-1))*(col+1);
    nnx = [pos I]';
    nnpr=net3(nnx);    
    sy = [xx;zz;nnpr];    
    vx=var(xx);
    vz=var(zz);
    mv=max(vx,vz);
    R=diag(repelem(mv+10e-6,2*ns+2));
    mkm = A*mk;
    Pkm = A*Pk*A' + Q;
    v = sy - H*mkm;
    S = H*Pkm*H' + R;
    K = (S'\(Pkm*H')')';
    mk = mkm + K*v;    
    Pk = Pkm - K*S*K';
    latcoords(k) = mk(1);  
    axcoords(k) = mk(2) - max(0,mk(3));
    elevcoords(k) = mk(4); 
    pos = [latcoords(k) axcoords(k)];
    
    
    if showplot
        IMAGE_env = abs(hilbert(images(:,:,k)'))';
        IMAGE_log = single(20*log10(IMAGE_env'));
        subplot(1,2,1)
        imagesc(ximg*1000,zimg*1000,IMAGE_log);
        caxis(mc - [DynRan 0]);
        line(1000*xpos,1000*zpos,'linestyle','none','marker','.','color','c','linewidth',4);
        line(1000*hydrophone(1),1000*hydrophone(3),'linestyle','none','marker','o','color','m','linewidth',2);
        set(gca,'fontsize',25,'fontweight','demi');
        xlabel('Lateral position [mm]');
        ylabel('Axial depth [mm]');
        colormap hot
        hold on
        plot(latcoords(k),axcoords(k),'Color','green','MarkerSize', 20,'Marker','.')
        ylim([zimg(1)*1000 zimg(end)*1000])
        xlim([ximg(1)*1000 ximg(end)*1000])
        if mk(4) < 1
            leg = 'Offset < 1 mm';
            color = 'green';
        else
            leg = append('Offset = ', string(round(mk(4),1)), ' mm');
            color = 'red';
        end
        legend(leg,'TextColor',color,'FontSize',10);
        hold off
        subplot(1,2,2)
        IMAGE_env = abs(hilbert(d));
        IMAGE_log = single(20*log10(IMAGE_env'))';
        imagesc(xm,ym,IMAGE_log);
        caxis(mc2 - [DynRan 0]);
        set(gca,'fontsize',25,'fontweight','demi');
        ylabel('Time [\mus]');
        yticks([0 5 10 15])
        xlabel('Source');
        title('Measured time series')
        drawnow
    end
end