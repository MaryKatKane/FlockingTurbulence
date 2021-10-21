clear
path(path,'./libDNS/');
fftw('planner','patient')
%% Parameters
%Basic definitions
global P_Phys;
global P_num;
P_Phys = struct();
P_num = struct();

% Physical parameters
P_Phys.V = 5;   % swimming velocity
P_Phys.R = 1; % radius of interaction
P_Phys.taup = .01; % response time of the angle
P_Phys.L = 2*pi; % domain size
P_Phys.nu = 5e-13; % hyperviscosity
P_Phys.p = 4; % power of Laplacian
P_Phys.mu = 2; % hypofriction
P_Phys.q = 4; % inverse power of Laplacian
P_Phys.kinf = 2; % minimal forcing wavenumber
P_Phys.ksup = 3; % maximal forcing wavenumber
P_Phys.stdf = 1; % forcing amplitude

% Numerical paremeters
P_num.n=256;
P_num.dt=1e-3; % tme step
P_num.Tmax=10; % total time of integration
P_num.Tsp=1e-2; % Spectral output every s
P_num.Tviz=1e-2; % Visualization of density every nviz steps
P_num.np = 1000;
P_num.srcsz = 0.01;


%% Definition of wavevectors
X1d=(0:(P_num.n-1))*(2*pi/P_num.n);
[Y,X]=meshgrid(X1d,X1d);
fpos=(0:(P_num.n/2));
fneg=(-(P_num.n/2-1)):-1;
kt=single([fpos,fneg]);
P_num.kx=meshgrid(kt)';
P_num.ky=meshgrid(kt);
P_num.ks=P_num.kx.*P_num.kx+P_num.ky.*P_num.ky;
clear('fpos','fneg','kt')
IIkill=getForcingIndex(P_num.n/3,P_num.n,P_num.ks,P_num.n); % For dealiasing
Kill=ones(P_num.n,P_num.n,'single');
Kill(IIkill)=0.0;

%% Initial conditions
Theta = 2*pi*rand(P_num.np,1);
Xp = 2*pi*rand(P_num.np,1);
Yp = 2*pi*rand(P_num.np,1);
load('initial_cond.mat','Omhat');

%% Compute spectra
P_num.KK = cell(P_num.n/2,1);
A = ceil(sqrt(P_num.ks));
for idk=1:(P_num.n/2-1)
    P_num.KK{idk} = find(A==idk-1);
end
P_num.KK{P_num.n/2} = find(A>=P_num.n/2-1);


om_sq = 1;
Ek=GetSpectra(calc_psi(Omhat));
Ekin = sum(Ek);
kk=0:((P_num.n/2)-1);

%% Plots
scrsz = get(groot,'ScreenSize');
ScreenTop=scrsz(4);
ScreenLeft=scrsz(3);

close all;

width=600*.85;
height=500*.85;

f1=figure(1);
FigureSize=[ScreenLeft-1.125*width ScreenTop-1*height,width,height];
set(gcf,'Position',FigureSize,'PaperPositionMode', 'auto')
pom=pcolor(X,Y,ifft2(Omhat,'symmetric'));
shading flat
hold on
hq = quiver(Xp,Yp,cos(Theta),sin(Theta),'Color','k');
colorbar
caxis([-30 30])
set(gca,'XTick',[],'YTick',[])
axis equal
axis tight
axis xy
xlim([0 P_Phys.L])
ylim([0 P_Phys.L])

f2=figure(2);
FigureSize=[ScreenLeft-2.25*width ScreenTop-1*height,width,height];
set(gcf,'Position',FigureSize,'PaperPositionMode', 'auto')
pEk=loglog(kk,Ek(end,:),'o-');
xlim([1 P_num.n/3])

f3=figure(3);
FigureSize=[ScreenLeft-2.25*width ScreenTop-2.5*height,width,height];
set(gcf,'Position',FigureSize,'PaperPositionMode', 'auto')
tout = 0;
pEtot=plot(tout,om_sq,'bo-');


%% Time loop
Diss = P_Phys.nu*double(P_num.ks).^P_Phys.p + P_Phys.mu./double(P_num.ks).^P_Phys.q;
ExpNu=Kill.*single(exp(-P_num.dt*Diss));
ExpNuh=Kill.*single(exp(-0.5*P_num.dt*Diss));
clear('Diss');
IIf=find(P_num.ks>=P_Phys.kinf^2 & P_num.ks<=P_Phys.ksup^2);
Force=zeros(P_num.n,P_num.n,'single');
stdforce = single(sqrt(P_num.dt)*P_num.n^2*P_Phys.stdf);
tviz = P_num.Tviz;
tsp = P_num.Tsp;
P_num.xgrid = [X X(:,1); X(1,:)+P_Phys.L X(1,1)+P_Phys.L];
P_num.ygrid = [Y Y(:,1)+P_Phys.L; Y(1,:) Y(1,1)+P_Phys.L];

testmov = false;
if testmov
    Mov = [];
end
%%
for t=P_num.dt:P_num.dt:P_num.Tmax
    
    % Add forcing (shot noise)
    Force(IIf) = stdforce.*(randn(size(IIf))+1i*randn(size(IIf)));
    Force=SymOp(Force);
    Omhat(IIf)=Omhat(IIf)+Force(IIf);
    
    % Runge-Kutta 2nd order
    [K,U,V] = RHS_omeg(Omhat);
    [Kx,Ky,Kthet] = RHS_part(double(U),double(V),double(ifft2(Omhat,'symmetric')),Xp,Yp,Theta);
    Om = ExpNuh.*(Omhat+.5*P_num.dt*K);
    Xptmp = mod(Xp+.5*P_num.dt*Kx,P_Phys.L);
    Yptmp = mod(Yp+.5*P_num.dt*Ky,P_Phys.L);
    Thetatmp = Theta+.5*P_num.dt*Kthet;
    [K,U,V] = RHS_omeg(Om);
    [Kx,Ky,Kthet] = RHS_part(double(U),double(V),double(ifft2(Om,'symmetric')),Xptmp,Yptmp,Thetatmp);
    Omhat=ExpNu.*(Omhat+P_num.dt*K);
    Xp = mod(Xp+P_num.dt*Kx,P_Phys.L);
    Yp = mod(Yp+P_num.dt*Ky,P_Phys.L);
    Theta = Theta+P_num.dt*Kthet;
    
    if isnan(Omhat(1,2))
        fprintf('Error: Theta is NaN\n');
        break;
    end
    
    % Visualization
    if t>=tviz
        set(hq,'XData',Xp,'YData',Yp,'UData',cos(Theta),'VData',sin(Theta));
        set(pom,'CData',ifft2(Omhat,'symmetric'))
        if testmov
            Mov = [Mov; getframe(1)];
        end
        pause(0.001)
        tviz = tviz + P_num.Tviz;
    end
    if t>=tsp
        Ekt=GetSpectra(calc_psi(Omhat));
        Ek=[Ek;Ekt];
        Ekin=[Ekin,sum(Ekt)];
        tout = [tout t];
        set(pEk,'YData',Ek(end,:));
        %set(pEtot,'XData',tout,'YData',Ekin);
        om_sq = [om_sq; mean(abs(Omhat.^2),'all')^0.5];
        set(pEtot,'Xdata',tout,'YData',om_sq);
        tsp=tsp+P_num.Tsp;    
    end
end
%%
if testmov
    myVideo = VideoWriter('flocks2d','MPEG-4');
    open(myVideo);
    writeVideo(myVideo,Mov);
    close(myVideo)
end