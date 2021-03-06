% file "plot_powspe.m"

clear all
close all

taum=0.01;
maxi=1000;

filename='data/2013-03-22_14-48_id-0_r-0.01_mu-1__5.dat';
%filename='data/2013-03-22_14-48_id-1_r-0.01_e-1__5.dat';
%filename='data/2013-03-22_14-48_id-2_r-0.01_cv-1__5.dat';
temp=importdata(filename,'\t',1);

tau=temp.data(1,5)*taum;
dt=temp.data(1,1)*taum;
N=temp.data(1,2);
df=1.0/(N*dt*taum); % [Hz]
fmax=1.0/(2.0*dt*taum); %[Hz]
f=[df:df:fmax]; % [Hz]
w=f*2.0*pi(); % [rad/s]
t=[0:dt:(dt*(N-1)/2)];

temp=importdata(filename,'',21);
S=temp.data();
C=df/2/pi*fft(S);
% plotting original data
figure(1)
loglog(f,S,'b-')
absC=abs(C);
figure(2)
semilogx(t(1:maxi),absC(1:maxi),'k-')
hold on
semilogx(t(1:maxi),real(C(1:maxi)),'b-')
semilogx(t(1:maxi),imag(C(1:maxi)),'r-')
hold off

% computing averaged data for given binsize
binsize=10;
f2=[(binsize/2*df):(binsize*df):(fmax-binsize/2*df)];
S2=zeros(1,length(f2));
std_dev_S2=zeros(1,length(f2));

for i=0:(length(f2)-1)
	tmp=S((i*binsize+1):(i+1)*binsize);
	S2(i+1)=mean(tmp);%/(1+4*pi^2*f2(i+1)^2*tau^2);
	std_dev_S2(i+1)=std(tmp);
end

% plotting averaged data
figure(3)
plot(f2(1:1000),S2(1:1000),'k-')

