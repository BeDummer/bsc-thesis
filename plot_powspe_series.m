% file "plot_powspe_series.m"

clear all

Gen=6;

filename='data/2013-02-19_14-40__';
filename_tmp=[filename '0.dat'];
temp=importdata(filename_tmp,'\t',1);

lines={'k-';'b-';'c-';'m-';'r-';'g-';'y-';'k--';'b--';'c--';'m--';'r--';'g--';'y--';'k:';'b:';'c:';'m:';'r:';'g:';'y:';'k-.';'b-.';'c-.';'m-.';'r-.';'g-.';'y-.'};
dt=temp.data(1,1);
N=temp.data(1,2);
df=1.0/(N*dt); % [Hz]
fmax=1.0/(2.0*dt); %[Hz]

binsize=10;
f2=[(binsize/2*df):(binsize*df):(fmax-binsize/2*df)];
S2=zeros(1,length(f2));
std_dev_S2=zeros(1,length(f2));

figure

for i=1:1:(Gen-1)
	filename_tmp=[filename int2str(i) '.dat'];
	temp=importdata(filename_tmp,'',7);
	S=temp.data();

	% computing averaged data for given binsize
	for j=0:(length(f2)-1)
		tmp=S((j*binsize+1):(j+1)*binsize);
		S2(j+1)=mean(tmp);
		std_dev_S2(j+1)=std(tmp);
	end
	% plotting averaged data
	plot(f2(1:200),S2(1:200),lines{i+1})
	if i==1
		hold on
	end

end
hold off
