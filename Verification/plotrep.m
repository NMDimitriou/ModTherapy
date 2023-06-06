function hh=plotrep(expdat,exp_ref,simfile,simdat,reps)

expd = readmatrix(['../Data/' expdat '.txt']);
expdr= readmatrix(['../Data/' exp_ref '.txt']);

inds = ismember(expdr(:,1),expd(:,1));

for i=1:reps
    simd{i} = readmatrix([simfile '/' simdat num2str(i) '.txt']);
end

run plotopt.m
f=figure('Position', [100, 100, 800, 700],'Visible','off');
hold on
for i=1:reps
    h{i}=plot(simd{i}(:,1),simd{i}(:,4),'b');
    h{i}.Color(4) = 0.01;  
end
errorbar(expdr(inds(:,1)==1,1),expdr(inds(:,1)==1,2),expdr(inds(:,1)==1,3),'.k','MarkerSize',20,'LineWidth',2)
errorbar(expdr(inds(:,1)==0,1),expdr(inds(:,1)==0,2),expdr(inds(:,1)==0,3),'.y','MarkerSize',20,'LineWidth',2)
xlabel('time (hours)')
ylabel('cell count')
axis([0 72 0 5e5])
%axis tight
hold off
print(f,['plot_' simfile '.png'],'-r350','-dpng')

end
