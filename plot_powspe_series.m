% file "plot_powspe_series.m"

clear all

start=0;
Gen=6;

filename='data/2013-03-08_16-09_0.3_0.4__';
filename_tmp=[filename '0.dat'];
temp=importdata(filename_tmp,'\t',1);

lines={'k-';'b-';'c-';'m-';'r-';'g-';'y-';'k--';'b--';'c--';'m--';'r--';'g--';'y--';'k:';'b:';'c:';'m:';'r:';'g:';'y:';'k-.';'b-.';'c-.';'m-.';'r-.';'g-.';'y-.'};
dt=temp.data(1,1);
N=temp.data(1,2);
k_max=temp.data(1,7)-1;
k=[1:k_max];
df=1.0/(N*dt); % [Hz]
fmax=1.0/(2.0*dt); %[Hz]

binsize=10;
f2=[(binsize/2*df):(binsize*df):(fmax-binsize/2*df)];
S2=zeros(1,length(f2));
std_dev_S2=zeros(1,length(f2));
rho_kum=zeros(1,Gen-1);
u_rho_kum=zeros(1,Gen-1);

for i=start:(Gen-1)
	filename_tmp=[filename int2str(i) '.dat'];

	temp=importdata(filename_tmp,'\t',4);
	r0=temp.data(1,2);
	
%	if i~=0
%		CV=temp.data(1,4);
%		u_CV=temp.data(2,4);
%	end

	temp=importdata(filename_tmp,'\t',9);
	rho=temp.data(1:k_max,1);
	std_dev=sqrt(temp.data(1:k_max,2));

	temp=importdata(filename_tmp,'',21);
	S=temp.data()./r0;

	% computing averaged data for given binsize
	for j=0:(length(f2)-1)
		tmp=S((j*binsize+1):(j+1)*binsize);
		S2(j+1)=mean(tmp);
		std_dev_S2(j+1)=std(tmp);
	end
	% plotting averaged data
	figure(1)
	loglog(f2(1:500),S2(1:500),lines{i+1})
	if i==start
		hold on
	end

	figure(2)
	errorbar(k,rho,std_dev,lines{i+1})
	if i==start
		hold on
	end

%	if i~=0
%		y0=mean(S(1:100));
%		u_y0=std(S(1:100));
%		rho_kum(i)=(y0/CV^2-1)/2;
%		u_rho_kum(i)=sqrt((u_y0/2/CV^2)^2+(u_CV*y0/CV^3)^2);
%	end
end

figure(1)
hold off
figure(2)
hold off

%rho_kum
%u_rho_kum
