%%
clc; clear
digits(16)

A=readmatrix('final_calibration_C33A_control_R1_werr_02.txt');
B=readmatrix('hist_calibration_C33A_control_R1_werr_02.txt');
C=readmatrix('rnd_calibration_C33A_control_R1_werr_02.txt');


figure;
hold on
plot(B(:,1),B(:,2),'LineWidth',2,'DisplayName','c++');
%histogram(A(:,1),256,'Normalization','probability','EdgeColor','none')
%[h,x]=ksdensity(A(:,1),'NumPoints',256);
%sh=sum(h);
%h=h/sh;
%plot(x,h,'LineWidth',2)

tmp = A(:,1);
dx=(max(tmp)-min(tmp))/256;
x_h = (min(tmp):dx:max(tmp));
p_hist = zeros(256,1);
for i=1:1024
    for k=1:256
        if(tmp(i)<=x_h(k+1) && tmp(i)>=x_h(k)) 
            p_hist(k) = p_hist(k) +1;
        end
    end
    %if(tmp(i)<=x_h(k+1) && tmp(i)>=x_h(k)) 
    %    p_hist(k) = p_hist(k) +1;
    %end
end

%normalize
p_hist = p_hist./sum(p_hist);

[N,edges,bin] = histcounts(A(:,1),'NumBins',256, 'BinLimits', [min(A(:,1)) max(A(:,1))], 'Normalization', 'probability');
plot(edges(1:end-1),N,'--','LineWidth',2,'DisplayName','matlab');
[hr,xr]=ksdensity(C(:,1),'NumPoints',256);
sr=sum(hr);
hr=hr/sr;
plot(xr,hr,'LineWidth',2,'DisplayName','ITS')
legend
axis tight
hold off

rel_err1 = sum(abs(N-B(:,2)'))/sum(B(:,2))

rel_err2 = sum(abs(p_hist-B(:,2)))/sum(B(:,2))

rel_err3 = sum(abs(p_hist-N'))/sum(N)