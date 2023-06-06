function hh=qoi_comp(expdat,path,cfile,vfile,sdat,reps)

expd = readmatrix(['../Data/' expdat '.txt']);

%inds = ismember(expdr,expd);

caldat=[];
valdat=[];

for i=1:reps
    tmp1   = readmatrix([path '/' cfile '/' sdat num2str(i) '.txt']);
    caldat = [caldat; tmp1(end,2)];
    tmp1   = readmatrix([path '/' vfile '/' sdat num2str(i) '.txt']);
    valdat = [valdat; tmp1(end,2)];
end
clear tmp1; 

run plotopt.m
disp('Plotting QOI posteriors');
f=figure('Position', [100, 100, 800, 700]); %,'Visible','off'
hold on
histogram(caldat,100,'EdgeColor','none','FaceAlpha',0.5,'Normalization','probability','DisplayName','Calibration dataset')
histogram(valdat,100,'EdgeColor','none','FaceAlpha',0.5,'Normalization','probability','DisplayName','Validation dataset')
stem(expd(end,2),.005,'k','LineWidth',2,'DisplayName','Observed')
xlabel('QOI')
ylabel('PDF')
%axis([-inf inf 0 1e6])
axis tight
legend('FontSize',20);
hold off
print(f,['plot_qoi_comp_' vfile '.png'],'-r350','-dpng')


disp('Performing KS test')
% cut off tails of predicted QOI probability distributions then
% compare using K-S test
% adopted from A. J. Connor. Calibration, Validation and Uncertainty Quantification. 2016.
m_cd=mean(caldat);
m_vd=mean(valdat);
std_cd=std(caldat);
std_vd=std(valdat);

pred_cal=caldat;
pred_val=valdat;

pred_cal(pred_cal < ( m_cd - 2*std_cd ))=[];
pred_cal(pred_cal > ( m_cd + 2*std_cd ))=[];
pred_val(pred_val < ( m_vd - 2*std_vd ))=[];
pred_val(pred_val > ( m_vd + 2*std_vd ))=[];
tol=0.05;
[h,p,ks2stat]=kstest2(pred_cal,pred_val,'Alpha',tol);

if(ks2stat < tol)
    A = ['K-S test accepted with ks2stat = ' num2str(ks2stat) ' < ' num2str(tol)];
else
    A = ['K-S test rejected with ks2stat = ' num2str(ks2stat) ' >= ' num2str(tol)];
end

% does predicted QOI agree with the observed ?

expSD  = expd(end,3);
diffSD = sqrt(std_vd.^2 + expSD.^2);

pDiscr = expd(end,2) - m_vd;

B = ['Error between predicted and observed QOI = ' num2str(abs(pDiscr)) ' \pm ' num2str(diffSD)];

writetable(cell2table({A;B}),['ks_qoi_' vfile '.txt']);
end

