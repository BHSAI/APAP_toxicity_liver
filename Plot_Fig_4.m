% Script to generate figure 4
% MFA data 10h acetaminophen
% Data is obtained from metabolic flux analysis in the current study
FN = {'G6PC','GPI','PYGL','ALDO','GAPDH^#','GK^#','ENO^#','PK','LDH','PC','PCK','CS','IDH','OGDH','PCC','SDH'};
% mean and std deviation for control samples
cmean = [49.26	1.35
38.64	2.35
10.62	1.62
38.64	2.35
30.82	1.69
7.82	1.07
30.82	1.69
30.08	6.07
46.55	2.33
76.63	7.88
91.71	9.07
43.56	3.94
43.56	3.94
43.56	3.94
15.08	1.55
58.64	5.41];
% APAP treatment values
tmean = [46.67	1.15
42.47	1.55
4.19	0.77
42.47	1.55
32.28	0.99
10.19	0.87
32.28	0.99
50.73	3.59
48.35	1.52
99.09	4.19
115.29	4.62
45.16	2.51
45.16	2.51
45.16	2.51
16.20	0.69
61.36	3.09];
mm = [cmean(:,1) tmean(:,1)]';
figure (4); set(figure(4),'Units','inches','Position',[0.5 0.5 8.8 7.0]);
subplot(2,1,1); subplot('position',[0.08 0.6 0.85 0.35]);
h = bar(mm','LineWidth',1.25); hold on;
set(h,'BarWidth',0.8); 
h(1).FaceColor = 'black';h(2).FaceColor = 'none';
N = nan(1,18);N(4)=18;N(9)=62;N(11)=110;N(12)=128;
plot(0:17,N,'k*','LineWidth',1.25,'MarkerSize',7)
ci = [cmean(:,2) tmean(:,2)]';
groupwidth = min(0.8, 2/(2+1.5));
for i = 1:2
x = (1:16) - groupwidth/2 + (2*i-1)*groupwidth/(2*2);
eh=errorbar(x,mm(i,:),ci(i,:),'k', 'linestyle', 'none','LineWidth',1.25);
end
xticks([1:16]); axis([0 17 0 140]);yticks([0:40:140]);
set(gca,'xticklabel',FN,'XTickLabelRotation',45)
ylabel('Flux (\mumol/kg/min)');
lg=legend(h,'Control (10h)','APAP (10h)'); legend('boxoff');
lg.FontSize = 11;
set(gca,'LineWidth',1.25,'FontName','Arial','FontSize',12);
text(-1.3,142,'a','FontName','Arial','FontSize',14)
% Aprroximate metabolic flux values obtained from the literature data for
% short term fasting conditions (5 h) 
% detailed in the main text of the manuscript
mean_5h = [59.11
29.56
29.56
29.56
39.41
19.70
39.41
19.23
29.76
49.00
58.64
27.85
27.85
27.85
9.64
37.5];
subplot(2,1,2); subplot('position',[0.08 0.11 0.85 0.35]);
h = bar(mean_5h,'LineWidth',1.25); hold on;
set(h,'BarWidth',0.25); 
h(1).FaceColor = 'black';xticks([1:16]);axis([0 17 0 140])
set(gca,'xticklabel',FN,'XTickLabelRotation',45)
ylabel('Flux (\mumol/kg/min)');
set(gca,'LineWidth',1.25,'FontName','Arial','FontSize',12);
yticks([0:40:140]);
lg=legend(h,'Control (5h)'); legend('boxoff');
lg.FontSize = 11;
text(-1.3,142,'b','FontName','Arial','FontSize',14)
