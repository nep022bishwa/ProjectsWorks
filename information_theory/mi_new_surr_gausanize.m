clear all; close all; clc;



clear;close all;
clc;
%reading of data 
load 'C:\Users\bishwa\Documents\Saturn_simulation_data\flow.mat';
load 'C:\Users\bishwa\Documents\Saturn_simulation_data\magnetic_field.mat';
load 'C:\Users\bishwa\Documents\Saturn_simulation_data\density.mat';
load 'C:\Users\bishwa\Documents\Saturn_simulation_data\co-ordinates.mat';

xe=xe(1:end-1,1:end-1);ye=-1*ye(1:end-1,1:end-1);


xe=double(xe);ye=double(ye);

k=700;
ne=den;
[NI,NJ,NT] = size(ne);

for i=1:NT
    V_y(:,:,i)=-1*V_y(:,:,i);
end

load '..\curl_without_background_flow_interp_cordinate.mat';


Rq=sqrt(xq.^2+yq.^2);
I_cond=(Rq<=7);

for i=1:NT
    my_val=curlz(:,:,i);
    
    my_val(I_cond)=nan;
    curlz(:,:,i)=my_val;
end


T0=1;                      
SimTime=(0:NT-1)*150;
figure()
pcolor(xe,ye,log10(ne(:,:,k)));
xlabel('R_s');ylabel('R_s');
shading interp
caxis([-3 0.5]);xlim([-100 30]),ylim([-60 60]);
daspect([1 1 1]);
h=colorbar;
ylabel(h, 'log_{10}(n_e)'); set(gca,'fontsize',30,'fontweight','bold');
set(gcf,'color','w'); hold on;
[nx,ny]=size(xq);

%%%%%% dawn longer box box 







%row=137; col=426; %This is closer box 
row=138; col=436;
ref_xx=xe(row:row+70,col:col+5);
ref_yy=ye(row:row+70,col:col+5);
top_x=[ref_xx(1,1),ref_xx(1,end),ref_xx(end,end),ref_xx(end,1),ref_xx(1,1)];
top_y=[ref_yy(1,1),ref_yy(1,end),ref_yy(end,end),ref_yy(end,1),ref_yy(1,1)];

plot(top_x,top_y,'k')




%getting the vorticity at the dusk region 

[in,on]=inpolygon(xq,yq,top_x,top_y);

 
  for tt=1:NT
    C=curlz(:,:,tt);
    top_BZ1= C(in); %electron density
      C1=randperm(length(top_BZ1));
  curl_top(:,tt)=(top_BZ1(C1(1:100)));
  
  end

  clear in on row col

  load '..\grid_cordinate.mat'


  plot(xx,yy,'k')



  %getting the vorticity at the grid box
for i=1:(length(xx))
 [in,on]=inpolygon(xq,yq,xx(:,i),yy(:,i));

 
  for tt=1:NT
    A=curlz(:,:,tt);
    BZ1=   A(in); %vorticity
     A1=randperm(length(BZ1));
  curl_grid(:,tt,i)=(BZ1(A1(1:100)));
 
  end
 
end



time=SimTime./3600;
II=find(time>10 & time<=time(end)-10);
%end_box_no=2;


nn=floor(length(II)./2);


ref_x=curl_top(:,1:nn);
ref_x=ref_x';

[ref_x] = gaussianize(ref_x); % This is gaussinize the data 
figure()
plot(ref_x,'.')

x=ref_x(:);
n_surrogates=50; 
%xsurr = multivariate_surrogates((curl_top(:,1:length(II)))',n_surrogates);

no_bin=2*ceil(log2(length(x))+1);
%no_bin=2*iqr(x )./nthroot(length(x),3);
 % number of surrogate data

box_no=5;   
curl_grid1=curl_grid(:,II,box_no);
curl_grid1=curl_grid1';

%ysurr = multivariate_surrogates(curl_grid1,n_surrogates);

%getting the mutual information 
for tt=1:nn
ref_y=curl_grid1(tt:nn+tt-1,:);
 ref_y = gaussianize(ref_y); % This also gausanize the data 
y=ref_y(:);

mi_value(tt)=mutual_information(x,y,no_bin);
   corr_mat=corr(ref_x,ref_y);
   lin_cor(tt)=trace(corr_mat)/100;
  % lin_cor(tt)=corr(ref_x,ref_y);

  % ygs = caaft(y, M,12*3600/150);

xsurr = multivariate_surrogates_new(ref_x,n_surrogates);
ysurr = multivariate_surrogates_new(ref_y,n_surrogates);
for i=1:n_surrogates
  ref_sur1=xsurr(:,:,i);
  ref_sur_x=ref_sur1(:);
  ref_sur2=ysurr(:,:,i);
   ref_sur_y=ref_sur2(:);
  mi_sur(i)=mutual_information(ref_sur_x,ref_sur_y,no_bin);
end
   
  mean_sur_mi(tt)=mean(mi_sur);
std_sur_mi(tt)=std(mi_sur);

end

upper_lim=mean_sur_mi+3*std_sur_mi;
lower_lim=mean_sur_mi-3*std_sur_mi;

figure()
subplot(2,1,1)
plot(time(1:nn)-time(1),mi_value,'b')
hold on 

plot(time(1:nn)-time(1), mean_sur_mi,'r')

plot(time(1:nn)-time(1), upper_lim,'g')
plot(time(1:nn)-time(1),lower_lim,'g');
%legend('MI','<MI_{surr}>','<MI_{surr}>+3\sigma','<MI_{surr}>-3\sigma')
set(gca,'fontsize',25);ylabel('MI')
title(strcat('Curlz vs curlz for box number=',num2str(box_no)))
hold off
subplot(2,1,2)
 plot(time(1:nn)-time(1),lin_cor,'b')
 set(gcf,'color','white'); xlabel('Time lag (hr)'); ylabel('Corr')
set(gca,'fontsize',25);legend('Corr')
set(gcf, 'Position', get(0, 'Screensize'));

big_lambda=sqrt(1-exp(-2*mi_value));


figure()

plot(time(1:nn)-time(1),mi_value,'b')
hold on 
plot(time(1:nn)-time(1),big_lambda,'r')
legend('I','\Lambda');set(gcf,'color','white');
set(gca,'fontsize',25);ylabel('MI')
title(strcat('Curlz vs curlz for box number=',num2str(box_no)))
%%%%%%%%%%%%% plot of surrogate 
figure()
subplot(2,2,1)
plot(time(1:nn),ref_x,'.b');hold on;plot(time(1:nn),xsurr(:,:,1),'.k');xlabel('Time lag (hr)'); 
subplot(2,2,2)
plot(time(1:nn),ref_x,'.b');hold on;plot(time(1:nn),xsurr(:,:,10),'.k');xlabel('Time lag (hr)'); 
subplot(2,2,3)
plot(time(1:nn),ref_x,'.b');hold on;plot(time(1:nn),xsurr(:,:,25),'.k');xlabel('Time lag (hr)'); 
subplot(2,2,4)
plot(time(1:nn),ref_x,'.b');hold on;plot(time(1:nn),xsurr(:,:,50),'.k');xlabel('Time lag (hr)');
legend('Data','Sur')


figure()
subplot(2,2,1)
plot(time(1:nn),ref_y,'.b');hold on;plot(time(1:nn),ysurr(:,:,1),'.k');xlabel('Time lag (hr)'); ylabel('Vorticity')
subplot(2,2,2)
plot(time(1:nn),ref_y,'.b');hold on;plot(time(1:nn),ysurr(:,:,10),'.k');xlabel('Time lag (hr)'); ylabel('Vorticity')
subplot(2,2,3)
plot(time(1:nn),ref_y,'.b');hold on;plot(time(1:nn),ysurr(:,:,25),'.k');xlabel('Time lag (hr)'); ylabel('Vorticity')
subplot(2,2,4)
plot(time(1:nn),ref_y,'.b');hold on;plot(time(1:nn),ysurr(:,:,50),'.k');xlabel('Time lag (hr)'); ylabel('Vorticity')

%%%%% cross-correlation 
figure(2)
maxlag = 12*3600/150;
subplot(2,2,1)
[cy lags] = xcorr(ref_y,maxlag,'normalized');


cs = [];
for j=1:n_surrogates;
     [c lags] = xcorr(ysurr(:,:,j),maxlag,'normalized');
     cs(:,:,j) = c;
end

csm = mean(cs,3);
cstd = std(cs,0,3);

kk=[1,2,3,101,102,103,201,202,203];
for j=1:9
subplot(3,3,j)
plot(lags,cy(:,kk(j)),lags,csm(:,kk(j)));
hold on
h = fill([lags fliplr(lags)],[transpose(csm(:,kk(j))+cstd(:,kk(j))) fliplr(transpose(csm(:,kk(j))-cstd(:,kk(j))))],'g');
h.FaceAlpha = 0.1;
%h.EdgeColor = 'none';
set(gca,'xlim',[-maxlag maxlag])
hold off
end
