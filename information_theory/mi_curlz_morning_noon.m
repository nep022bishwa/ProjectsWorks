clear all;
clc; 
close all;


clear;close all;
clc;
%reading of data 
load 'C:\Users\bishwa\Documents\Saturn_simulation_data\flow.mat';
load 'C:\Users\bishwa\Documents\Saturn_simulation_data\magnetic_field.mat';
load 'C:\Users\bishwa\Documents\Saturn_simulation_data\density.mat';
load 'C:\Users\bishwa\Documents\Saturn_simulation_data\co-ordinates.mat';

xe=xe(1:end-1,1:end-1);ye=-1*ye(1:end-1,1:end-1);


xe=double(xe);ye=double(ye);

k=900;
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
pcolor(xq,yq,curlz(:,:,k));
xlabel('R_s');ylabel('R_s');
shading interp
caxis([-150 150]);xlim([-100 30]),ylim([-60 60]);
daspect([1 1 1]);
h=colorbar;
ylabel(h, 'Curl V'); set(gca,'fontsize',30,'fontweight','bold');
set(gcf,'color','w'); hold on;
[nx,ny]=size(xq);


%%%%%% dawn longer box box 
%row=137; col=426; %This is closer box 
for i=1:NI
    for j=1:NJ
theta(i,j)=angle(xe(i,j),ye(i,j));
    end
end
R=sqrt(xe.^2+ye.^2);


%getting the box for dawn region
[row1,col1]=find((R>22.0 & R<22.50) & (theta>250 & theta<250.5));
%row=137; col=426; %This is closer box 
row=row1(1); col=col1(1);



ref_xx=xe(row:row+10,col:col+30);
ref_yy=ye(row:row+10,col:col+30);
dawn_x=[ref_xx(1,1),ref_xx(1,end),ref_xx(end,end),ref_xx(end,1),ref_xx(1,1)];
dawn_y=[ref_yy(1,1),ref_yy(1,end),ref_yy(end,end),ref_yy(end,1),ref_yy(1,1)];

plot(dawn_x,dawn_y,'k')

clear in on row col

  load '..\grid_cordinate.mat'


  plot(xx,yy,'k')





  %getting the vorticity at the dawn region 

[in,on]=inpolygon(xq,yq,dawn_x,dawn_y);

% no_points=400;
  for tt=1:NT
    C=curlz(:,:,tt);
    
    dawn_curlz= C(in); %electron density
    
     
    if tt==1;
        
      dawn_random=randperm(length(dawn_curlz));
     no_points=length(dawn_curlz);
    end
     
 % curl_top(tt,:)=(top_BZ1(C1(1:no_points)));
  %curl_dawn(tt,:)=dawn_curlz(dawn_random(1:no_points));
  curl_dawn(tt,:)=dawn_curlz(dawn_random(1:no_points));
  end


  clear in on 

  %getting the box for noon region 
[row1,col1]=find((R>17.5 & R<17.75) & (theta>30 & theta<30.5));
%row=137; col=426; %This is closer box 
row=row1(1); col=col1(1);
ref_xx=xe(row:row+25,col:col+50);
ref_yy=ye(row:row+25,col:col+50);
noon_x=[ref_xx(1,1),ref_xx(1,end),ref_xx(end,end),ref_xx(end,1),ref_xx(1,1)];
noon_y=[ref_yy(1,1),ref_yy(1,end),ref_yy(end,end),ref_yy(end,1),ref_yy(1,1)];

plot(noon_x,noon_y,'k')




%getting the vorticity at the dusk region 

[in,on]=inpolygon(xq,yq,noon_x,noon_y);

% no_points=400;
  for tt=1:NT
    C=curlz(:,:,tt);
    
    noon_curlz= C(in); %electron density
    
    if tt==1;
        %this gives the random permutation of the noon curlz
      noon_random=randperm(length(noon_curlz));
    end
    
 % curl_top(tt,:)=(top_BZ1(C1(1:no_points)));
  curl_noon(tt,:)=noon_curlz(noon_random(1:no_points));
  end



  time=SimTime./3600;
II=find(time>=10 & time<=time(end));
%end_box_no=2;


nn=floor(length(II)./2);


curl_dawn1=curl_dawn(II,:);

curl_dawn1=curl_dawn1-mean(curl_dawn1);

n_surrogates=200; % This gives the number of surrogate


%no_bin=ceil(log2(length(x))+1);

no_bin=30; % This number of bins is used for the variation of the
%vorticity MI calculation

%no_bin=10;

curl_noon1=curl_noon(II,:);

curl_noon1=curl_noon1-mean(curl_noon1);
%ysurr = multivariate_surrogates(curl_grid1,n_surrogates);

yy=[curl_dawn1 curl_noon1];

figure()
subplot(2,1,1)
plot(time(II),curl_dawn1,'.');
xlabel('Time[Hr]')
subplot(2,1,2)
plot(time(II),curl_noon1,'.');
xlabel('Time [hr]')

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
title('Curlz vs curlz for dawn and noon box')
hold off
subplot(2,1,2)
 plot(time(1:nn)-time(1),lin_cor,'b'); hold on 
 plot(time(1:nn)-time(1), mean_sur_cor,'r')
 plot(time(1:nn)-time(1), upper_lim_cor,'g')
plot(time(1:nn)-time(1),lower_lim_cor,'g')
 set(gcf,'color','white'); xlabel('Time lag (hr)'); ylabel('Corr')
set(gca,'fontsize',25);legend('Corr')
set(gcf, 'Position', get(0, 'Screensize'));

%get the significane 

tot_signi=(mi_value-mean_sur_mi)./std_sur_mi;

figure()
plot(time(1:nn)-time(1),tot_signi,'linewidth',3)
set(gca,'fontsize',25);set(gcf,'color','w');
xlabel('Time[hr]');ylabel('Significnace')

%getting the analysis for max mi

I_max=102;
%I_max=find(abs(mi_value & time(II)<6)==max(abs(mi_value)));


 y_value=yy(I_max:end,no_points+1:end);
 x_value=yy(1:end-I_max+1,1:no_points);

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
xlabel('Vort'); ylabel('vort');
set(gcf,'Color','w');set(gca,'fontsize',15)
hold off
