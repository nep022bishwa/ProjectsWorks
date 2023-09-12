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
  curl_top(tt)=mean(top_BZ1(C1(1:no_points)));
  
  end


  time=SimTime./3600;
II=find(time>=5 & time<=time(end));
%end_box_no=2;


nn=floor(length(II)./2);


curl_top1=curl_top(II);

curl_top1=curl_top1-mean(curl_top1(:));

n_surrogates=50; % This gives the number of surrogate


%no_bin=ceil(log2(length(x))+1);

no_bin=10; % This number of bins is used for the variation of the


yy=[curl_top1' curl_top1'];

%no_bin=ceil(log2(length(yy))+1);

for tt=1:nn
    %{
   ref_x= yy(1:end-tt+1,1);
   x=ref_x(:);
   ref_y=yy(tt:end,2);
   y=ref_y(:);
    %}
    %%{
ref_x=yy(1:nn,1);
x=ref_x(:);
ref_y=yy(tt:nn+tt-1,2);
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
   ref_sur1= ysurr(1:end-tt+1,1,i);
  ref_sur_x=ref_sur1(:);
   ref_sur2=ysurr(tt:end,2,i);
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
set(gca,'YScale','log')
%legend('MI','<MI_{surr}>','<MI_{surr}>+3\sigma','<MI_{surr}>-3\sigma')
set(gca,'fontsize',25);ylabel('MI')
title('Self Mutual information of box on dusk side')
hold off
subplot(2,1,2)
 plot(time(1:nn)-time(1),lin_cor,'b'); hold on 
 plot(time(1:nn)-time(1), mean_sur_cor,'r')
 plot(time(1:nn)-time(1), upper_lim_cor,'g')
plot(time(1:nn)-time(1),lower_lim_cor,'g')
 set(gcf,'color','white'); xlabel('Time lag (hr)'); ylabel('Corr')
set(gca,'fontsize',25);legend('Corr')

 signi=(mi_value-mean_sur_mi)./std_sur_mi;

 figure()
plot(time(1:nn)-time(1),signi,'LineWidth',3)
set(gca,'fontsize',15,'fontweight','bold')
set(gcf,'Color','white');
xlabel('Time [hr]'); ylabel('Significane')
set(gca,'fontsize',25);