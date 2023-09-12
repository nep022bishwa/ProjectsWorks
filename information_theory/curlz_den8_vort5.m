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


clear in on row col

  load '..\grid_cordinate.mat'


  plot(xx,yy,'k')



 %getting the vorticity at the grid box
box_no1=5;

 [in,on]=inpolygon(xq,yq,xx(:,box_no1),yy(:,box_no1));

 
  for tt=1:NT
    A=B_z(:,:,tt);
    BZ1=   A(in); %vorticity
    I_neg=find(BZ1<0);
   no_events(tt,:)=length(I_neg)/length(BZ1);
  end
 
  %getting the density in the box 5
   box_no2=8;
  [in,on]=inpolygon(xq,yq,xx(:,box_no2),yy(:,box_no2));

 
  for tt=1:NT
    A1=ne(:,:,tt);
    NE2=   A1(in); %vorticity
    Den_box(tt,:)=std(NE2);
  end



time=SimTime/3600;
figure()
subplot(2,1,1)
plot(time,100*no_events)
xlabel('Time'); ylabel('% BZ<0 events')
subplot(2,1,2)
plot(time,Den_box)
xlabel('Time'); ylabel('Fluctuation of Density')