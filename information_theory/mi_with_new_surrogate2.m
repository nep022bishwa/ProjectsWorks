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

 no_points=400;
  for tt=1:NT
    C=curlz(:,:,tt);
    
    top_BZ1= C(in); %electron density
    if tt==1;
      C1=randperm(length(top_BZ1));
    end
   mean_val=mean(top_BZ1);
   %curl_top(tt,:)= sqrt((top_BZ1(C1(1:no_points))-mean_val).^2);
  curl_top(tt,:)=(top_BZ1(C1(1:no_points)));
  
  end

  clear in on row col

  load '..\grid_cordinate.mat'


  plot(xx,yy,'k')



  %getting the vorticity at the grid box
box_no=5;

 [in,on]=inpolygon(xq,yq,xx(:,box_no),yy(:,box_no));

 
  for tt=1:NT
    A=curlz(:,:,tt);
    BZ1=   A(in); %vorticity

    if tt==1;
      A1=randperm(length(BZ1));
    end
     
    
     mean_val=mean(BZ1(:));
   curl_grid(tt,:)=sqrt((BZ1(A1(1:no_points))-mean_val).^2);
 %curl_grid(tt,:)=BZ1(A1(1:no_points));
  end
 




time=SimTime./3600;
II=find(time>=10 & time<=time(end));
%end_box_no=2;


nn=floor(length(II)./2);


curl_top1=(curl_top(II,:));



n_surrogates=50; % This gives the number of surrogate


%no_bin=ceil(log2(length(x))+1);

no_bin=40; % This number of bins is used for the variation of the
%vorticity MI calculation

%no_bin=10;

curl_grid1=log(curl_grid(II,:));


%ysurr = multivariate_surrogates(curl_grid1,n_surrogates);

yy=[curl_top1 curl_grid1];


%yy=[curl_top1 curl_top1];


%ref_x=yy(1:nn,1:no_points);
%
% x=ref_x(:);



%getting the mutual information 
for tt=1:nn
    %%{
   ref_x= yy(1:end-tt+1,1:no_points);
   x=ref_x(:);
   ref_y=yy(tt:end,no_points+1:end);
   y=ref_y(:);
    %}
    %{
ref_y=yy(tt:nn+tt-1,no_points+1:end);
y=ref_y(:);
    %}
%no_bin=ceil(log2(length(x))+1);
mi_value(tt)=mutual_information(x,y,no_bin);
   lin_cor(tt)=corr(x,y);

%small_lambda(tt)= sqrt(1-det(cov(x,y))/(var(x).*var(y)));

%getting the surrogate data
  ysurr = multivariate_surrogates_fix(yy,n_surrogates);

for i=1:n_surrogates

%%{
   ref_sur1= ysurr(1:end-tt+1,1:no_points,i);
  ref_sur_x=ref_sur1(:);
   ref_sur2=ysurr(tt:end,no_points+1:end,i);
   ref_sur_y=ref_sur2(:);
    %}

%{
  ref_sur1=ysurr(1:nn,1:no_points,i);
  ref_sur_x=ref_sur1(:);
  
 % ref_sur2=ysurr(1:nn,no_points+1:end,i);

 ref_sur2=ysurr(tt:nn+tt-1,no_points+1:end,i);
 ref_sur_y=ref_sur2(:);
%}
  mi_sur(i)=mutual_information(ref_sur_x,ref_sur_y,no_bin);
  lin_cor_sur(i)=corr(ref_sur_x,ref_sur_y);
end
   
  mean_sur_mi(tt)=mean(mi_sur);
std_sur_mi(tt)=std(mi_sur);

mean_sur_cor(tt)=mean(lin_cor_sur);
std_sur_cor(tt)=std(lin_cor_sur);
end

big_lambda=sqrt(1-exp(-2*mi_value));

upper_lim=mean_sur_mi+3*std_sur_mi;
lower_lim=mean_sur_mi-3*std_sur_mi;


upper_lim_cor=mean_sur_cor+3*std_sur_cor;
lower_lim_cor=mean_sur_cor-3*std_sur_cor;

figure()
subplot(2,1,1)
plot(time(1:nn)-time(1),mi_value,'b')
hold on 

plot(time(1:nn)-time(1), mean_sur_mi,'r')

plot(time(1:nn)-time(1), upper_lim,'g')
plot(time(1:nn)-time(1),lower_lim,'g');
%set(gca,'yscale','log')
%legend('MI','<MI_{surr}>','<MI_{surr}>+3\sigma','<MI_{surr}>-3\sigma')
set(gca,'fontsize',25);ylabel('MI')
title(strcat('Curlz vs curlz for box number=',num2str(box_no)))
hold off
subplot(2,1,2)
 plot(time(1:nn)-time(1),lin_cor,'b'); hold on 
 plot(time(1:nn)-time(1), mean_sur_cor,'r')
 plot(time(1:nn)-time(1), upper_lim_cor,'g')
plot(time(1:nn)-time(1),lower_lim_cor,'g')
 set(gcf,'color','white'); xlabel('Time lag (hr)'); ylabel('Corr')
set(gca,'fontsize',25);legend('Corr')
set(gcf, 'Position', get(0, 'Screensize'));
%{
figure()
plot(time(1:nn)-time(1),small_lambda,'y');hold on
plot(time(1:nn)-time(1),lin_cor,'--b');
plot(time(1:nn)-time(1),big_lambda,'r')
 set(gcf,'color','white'); xlabel('Time lag (hr)'); 
set(gca,'fontsize',15);legend('\lambda','corr','\Lambda')
set(gcf, 'Position', get(0, 'Screensize'));
%}
%getting the significance value for the maximum mi value time series 

 I=find(abs(mi_value)==max(abs(mi_value)));


 y_value=yy(I:end,no_points+1:end);
 x_value=yy(1:end-I+1,1:no_points);

no_bin=5:100;

y_max=[x_value y_value];
ref_y_max=y_value(:);
ref_x_max=x_value(:);

ysur_max=multivariate_surrogates_fix(y_max,n_surrogates);

for i=1:length(no_bin)
    
mi_value_max(i)=mutual_information(ref_x_max,ref_y_max,no_bin(i));

for j=1:n_surrogates
    ref_sur_x_max=ysur_max(:,1:no_points,j);
ref_sur_y_max=ysur_max(:,no_points+1:end,j);
sur_mi_max(j)=mutual_information(ref_sur_x_max(:),ref_sur_y_max(:),no_bin(i));
end
mean_sur_max(i)=mean(sur_mi_max);
std_sur_max(i)=std(sur_mi_max);
end


figure()
plot(no_bin,mi_value_max)

sig=abs(mi_value_max-mean_sur_max)./std_sur_max;

figure()

plot(no_bin,sig)
xlabel('# bin'); ylabel('significance'); 
set(gcf, 'color','w')
hold off


figure()
plot(ref_x_max,ref_y_max,'.')
hold on;
haxes = dscatter(ref_x_max,ref_y_max,'PLOTTYPE','contour');
xlabel('Vort'); ylabel('log(var(vort))');
set(gcf,'Color','w');set(gca,'fontsize',15)
hold off
%set(gca,'YScale','log')
%Fitting of the normal variables [ref_x and ref_y];
figure()

X=[ref_x_max ref_y_max];
[sigma,mu] = robustcov(X);
%sigma=cov(X); mu=mean(X);
Y = mvnrnd(mu,sigma,length(ref_x_max));





xmax=max(max(ref_x_max),max(Y(:,1)));
xmin=min(min(ref_x_max),min(Y(:,1)));

XEDGES=linspace(xmin,xmax,50);

ymax=max(max(ref_y_max),max(Y(:,2)));
ymin=min(min(ref_y_max),min(Y(:,2)));

YEDGES=linspace(ymin,ymax,50);


h1=histogram2(ref_x_max,ref_y_max,XEDGES,YEDGES,'Normalization','pdf');
l1=h1.Values;

h2=histogram2(Y(:,1),Y(:,2),XEDGES,YEDGES,'Normalization','pdf');
l2=h2.Values;

val2=l1-l2;


figure()
subplot(1,3,1)
surf(l1)

subplot(1,3,2)
surf(l2)
subplot(1,3,3)
surf(val2)


%%% getting the indices for maximum std 
%{

I_max_std=find(upper_lim==max(upper_lim));


y_std_max_value=yy(I_max_std:nn+I_max_std-1,no_points+1:end);
y_std_max_value=y_std_max_value(:);
x_std_max_value=ref_x(:);

figure()

plot(x_std_max_value,y_std_max_value,'.')
hold on;
haxes = dscatter(x_std_max_value,y_std_max_value,'PLOTTYPE','contour');
xlabel('Vort'); ylabel('log(var(vort))');
set(gcf,'Color','w');set(gca,'fontsize',15);
title('distribution corresponding to maximum standard deviation');
hold off

%}





%{
%%%% plot of surrogate 
figure()
subplot(2,2,1)
plot(time(1:nn),ref_x,'.b');hold on;plot(time(1:nn),ysurr(:,1:100,1),'.k')
subplot(2,2,2)
plot(time(1:nn),ref_x,'.b');hold on;plot(time(1:nn),ysurr(:,1:100,10),'.k')
subplot(2,2,3)
plot(time(1:nn),ref_x,'.b');hold on;plot(time(1:nn),ysurr(:,1:100,25),'.k')
subplot(2,2,4)
plot(time(1:nn),ref_x,'.b');hold on;plot(time(1:nn),ysurr(:,1:100,50),'.k')


figure()
subplot(2,2,1)
plot(time(1:nn),ref_y,'.b');hold on;plot(time(1:nn),ysurr(:,101:200,1),'.k')
subplot(2,2,2)
plot(time(1:nn),ref_y,'.b');hold on;plot(time(1:nn),ysurr(:,101:200,10),'.k')
subplot(2,2,3)
plot(time(1:nn),ref_y,'.b');hold on;plot(time(1:nn),ysurr(:,101:200,25),'.k')
subplot(2,2,4)
plot(time(1:nn),ref_y,'.b');hold on;plot(time(1:nn),ysurr(:,101:200,50),'.k')

%%%%% cross-correlation 
figure()
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
figure()
for j=1:9
subplot(3,3,j)
plot(lags,cy(:,kk(j)),lags,csm(:,kk(j)))
hold on
h = fill([lags fliplr(lags)],[transpose(csm(:,kk(j))+cstd(:,kk(j))) fliplr(transpose(csm(:,kk(j))-cstd(:,kk(j))))],'g');
h.FaceAlpha = 0.1;
%h.EdgeColor = 'none';
set(gca,'xlim',[-maxlag maxlag])
hold off
end

figure()
subplot(2,1,1)
plot(time(1:nn),ref_x,'.b'); hold on;plot(time(1:nn),mean(ref_x,2),'y','LineWidth',5)
subplot(2,1,2)
plot(time(1:nn),xsurr(:,:,25),'.k');hold on;plot(time(1:nn),mean(xsurr(:,:,25),2),'y','LineWidth',5)
figure()
subplot(2,1,1)
plot(time(1:nn),ref_y,'.b'); hold on;plot(time(1:nn),mean(ref_y,2),'y','LineWidth',5)
subplot(2,1,2)
plot(time(1:nn),ysurr(:,:,25),'.k');hold on;plot(time(1:nn),mean(ysurr(:,:,25),2),'y','LineWidth',5)

%plot(ref_x,curl_grid1(188:nn+188-1),'.','markersize',10)
%}