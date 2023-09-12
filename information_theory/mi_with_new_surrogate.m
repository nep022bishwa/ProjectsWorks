clear all;
close all;
clc;

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
figure(2)
pcolor(xe,ye,log10(ne(:,:,k)));
xlabel('R_s');ylabel('R_s');
shading interp
caxis([-3 0.5]);xlim([-100 30]),ylim([-60 60]);
daspect([1 1 1]);
h=colorbar;
ylabel(h, 'log_{10}(n_e)'); set(gca,'fontsize',30,'fontweight','bold');
set(gcf,'color','w'); hold on;
[nx,ny]=size(xq);

%%%%%% dawn box box 
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
  curl_top(:,tt)=(top_BZ1(C1(1:50)));
  
  end

  clear in on row col

  load '..\grid_cordinate.mat'

%%%%%for the grid box 
%{
row=84; col=304;
ref_row=row(1);

% give me the row and colum indices for the box on starting region
nx_box=3; ny_box=3;
xx=[]; yy=[];
for ii=1:nx_box
    ref_col=col(1);
    for jj=1:ny_box
       
start_ind_xe(ii,jj)=[ref_row(1)];
start_ind_ye(ii,jj)=[ref_col(1)];
ref_col=ref_col+40;
    end
    ref_row=ref_row+40;
end
start_ind_xe1=start_ind_xe(:);
start_ind_ye1=start_ind_ye(:);


for i=1:(nx_box*ny_box)
ref_xx=xe(start_ind_xe1(i):start_ind_xe1(i)+40,start_ind_ye1(i):start_ind_ye1(i)+40);
xx(:,i)=[ref_xx(1,1),ref_xx(1,end),ref_xx(end,end),ref_xx(end,1),ref_xx(1,1)];
ref_yy=ye(start_ind_xe1(i):start_ind_xe1(i)+40,start_ind_ye1(i):start_ind_ye1(i)+40);
yy(:,i)=[ref_yy(1,1),ref_yy(1,end),ref_yy(end,end),ref_yy(end,1),ref_yy(1,1)];
end
%}
plot(xx,yy,'k')



%getting the vorticity at the grid box
for i=1:(length(xx))
 [in,on]=inpolygon(xq,yq,xx(:,i),yy(:,i));

 
  for tt=1:NT
    A=curlz(:,:,tt);
    BZ1=   A(in); %vorticity
     A1=randperm(length(BZ1));
  curl_grid(:,tt,i)=(BZ1(A1(1:50)));
 
  end
 
end




time=SimTime./3600;
II=find(time>10 & time<=time(end)-10);
%end_box_no=2;


nn=floor(length(II)./2);


ref_x=curl_top(:,1:nn);

x=ref_x(:);
%x=sign(x).*log2(abs(x));
%curl_grid1=curl_grid(:,II,:);

no_bin=ceil(log2(length(x))+1);

%this is the calculation of mutal information 
for i=5%1:length(xx)
    
curl_grid1=curl_grid(:,II,i);

%getting the mutual information 
for tt=1:nn
ref_y=curl_grid1(:,tt:nn+tt-1);
y=ref_y(:);
%y=sign(y).*log2(abs(y));
mi_value(tt,i)=mutual_information(x,y,no_bin);
   lin_cor(tt,i)=corr(x,(y));

for j=1:50
      
   
       ygs = get_surrogate(y');
      
       sur_mi(j)=mutual_information(x,ygs,no_bin);
   end
  mean_sur_mi(tt,i)=mean(sur_mi);
std_sur_mi(tt,i)=std(sur_mi);

end
i

end





for i=5%1:9
    figure()
subplot(2,1,1)
plot(time(1:nn)-time(1),mi_value(:,i),'b')
hold on 

plot(time(1:nn)-time(1), mean_sur_mi(:,i),'r')

plot(time(1:nn)-time(1), mean_sur_mi(:,i)+3*std_sur_mi(:,i),'g')
plot(time(1:nn)-time(1), mean_sur_mi(:,i)-3*std_sur_mi(:,i),'g')
%legend('MI','<MI_{surr}>','<MI_{surr}>+3\sigma','<MI_{surr}>-3\sigma')
set(gca,'fontsize',25);ylabel('MI')
title(strcat('Curlz vs curlz for box number=',num2str(i)))
hold off
subplot(2,1,2)
 plot(time(1:nn)-time(1),lin_cor(:,i),'b')
 set(gcf,'color','white'); xlabel('Time lag (hr)'); ylabel('Corr')
set(gca,'fontsize',25);legend('Corr')
set(gcf, 'Position', get(0, 'Screensize'));
%saveas(gcf,strcat(num2str(i),'curlz_curlz_withour_background_flow.jpg'));
end

%xedge=logspace(-5,3,N); xedge=sort([-xedge xedge]);