%% Figure 2a, b: plot differential gene expression against q values: volcano plots
% load the differential gene expression data along with q values
% Tables APAP5h.txt and APAP10h.txt were copied from the
% Sleuth_table_liver files
% Remove "Insertpath" portion and include the path to the folder Weights_apap located on your machine
Gene_5h=readtable('"Insertpath"\Weights_apap\APAP5h.txt','Delimiter','\t');
Gene_10h=readtable('"Insertpath"\Weights_apap\APAP10h.txt','Delimiter','\t');
G5h_x = exp(Gene_5h{:,2});G5h_x = log2(G5h_x);
G5h_y = Gene_5h{:,1}; 
G10h_x = exp(Gene_10h{:,2});G10h_x = log2(G10h_x);
G10h_y = Gene_10h{:,1};
G10h_y(1) = 8.61e-55;
figure (1); set(figure(1),'Units','inches','Position',[0.5 0.5 9.8 7.0]);
subplot(2,2,1); subplot('position',[0.08 0.58 0.4 0.38]);
hold on; box off; grid off;
scatter(G5h_x(G5h_y<0.1 & G5h_x>0),-log10(G5h_y(G5h_y<0.1 & G5h_x>0)),20,'ro','filled');
scatter(G5h_x(G5h_y<0.1 & G5h_x<0),-log10(G5h_y(G5h_y<0.1 & G5h_x<0)),20,'go','filled');
scatter(G5h_x(G5h_y>=0.1),-log10(G5h_y(G5h_y>=0.1)),20,'ko','filled')
set(gca,'LineWidth',1.25,'FontName','Arial','FontSize',12);
xlabel('log_2(Fold change)'); ylabel('-log_1_0(FDR)');
text(-7.9,62,'a','FontName','Arial','FontSize',14)
axis([-6 6 0 60]); xticks([-6:3:6]);
yticklabels = get(gca, 'YTickLabel');
yticklabels{2}='';yticklabels{4}='';yticklabels{6}='';
set(gca, 'YTickLabel', yticklabels);
lg=text(3.35,62,'APAP (5 h)');
lg.FontSize = 11;lg.FontName = 'arial';

subplot(2,2,2); subplot('position',[0.57 0.58 0.4 0.38]);
hold on; box off; grid off;
scatter(G10h_x(G10h_y<0.1 & G10h_x>0),-log10(G10h_y(G10h_y<0.1 & G10h_x>0)),20,'ro','filled');
scatter(G10h_x(G10h_y<0.1 & G10h_x<0),-log10(G10h_y(G10h_y<0.1 & G10h_x<0)),20,'go','filled');
scatter(G10h_x(G10h_y>=0.1),-log10(G10h_y(G10h_y>=0.1)),20,'ko','filled')
% plot(G10h_x,-log10(G10h_y),'bo')
set(gca,'LineWidth',1.25,'FontName','Arial','FontSize',12);
xlabel('log_2(Fold change)'); ylabel('-log_1_0(FDR)');
text(-7.9,62,'b','FontName','Arial','FontSize',14)
axis([-6 6 0 60]); xticks([-6:3:6]);
yticklabels = get(gca, 'YTickLabel');
yticklabels{2}='';yticklabels{4}='';yticklabels{6}='';
set(gca, 'YTickLabel', yticklabels);

lg=text(3.1,62,'APAP (10 h)');
lg.FontSize = 11;lg.FontName = 'arial';

%% Figure 2c,d: plot differential gene expression against q values: volcano plots for genes mapped to iRno model
% load the differential gene expression data along with q values
% Tables for metabolic genes were obtained after mapping genes to the iRno model from Sleuth_table_liver files
% Remove "Insertpath" portion and include the path to the folder Weights_apap located on your machine
mGene_5h=tdfread('"Insertpath"\Weights_apap\Metabolic_genes_liver5h.tsv','\t');
mGene_10h=tdfread('"Insertpath"\Weights_apap\Metabolic_genes_liver10h.tsv','\t');
mG5h_x = exp(mGene_5h.b);mG5h_x = log2(mG5h_x);
G5h_y = mGene_5h.qval; 
mG10h_x = exp(mGene_10h.b);mG10h_x = log2(mG10h_x);
mG10h_y = mGene_10h.qval;
mG10h_y(251) = 8.61e-55; % oriignal value is 8.6081e-81, only reduced to fit in the axis scale
% figure (1); set(figure(1),'Units','inches','Position',[0.5 0.5 9.8 7.0]);
subplot(2,2,3); subplot('position',[0.08 0.1 0.4 0.38]);
hold on; box off; grid off;
scatter(mG5h_x(G5h_y<0.1 & mG5h_x>0),-log10(G5h_y(G5h_y<0.1 & mG5h_x>0)),20,'ro','filled');hold on;
scatter(mG5h_x(G5h_y<0.1 & mG5h_x<0),-log10(G5h_y(G5h_y<0.1 & mG5h_x<0)),20,'go','filled');hold on; box off
scatter(mG5h_x(G5h_y>=0.1),-log10(G5h_y(G5h_y>=0.1)),20,'ko','filled')
% plot(G5h_x,-log10(G5h_y),'ro')
set(gca,'LineWidth',1.25,'FontName','Arial','FontSize',12);
xlabel('log_2(Fold change)'); ylabel('-log_1_0(FDR)');
text(-7.9,62,'c','FontName','Arial','FontSize',14)
axis([-6 6 0 60])
xticks([-6:3:6]);
yticklabels = get(gca, 'YTickLabel');
yticklabels{2}='';yticklabels{4}='';yticklabels{6}='';
set(gca, 'YTickLabel', yticklabels);
% lg=text(3.35,56,'APAP (5h)');
lg.FontSize = 11;lg.FontName = 'arial';

subplot(2,2,4); subplot('position',[0.57 0.1 0.4 0.38]);
hold on; box off; grid off;
scatter(mG10h_x(mG10h_y<0.1 & mG10h_x>0),-log10(mG10h_y(mG10h_y<0.1 & mG10h_x>0)),20,'ro','filled');hold on; box off
scatter(mG10h_x(mG10h_y<0.1 & mG10h_x<0),-log10(mG10h_y(mG10h_y<0.1 & mG10h_x<0)),20,'go','filled');hold on; box off
scatter(mG10h_x(mG10h_y>=0.1),-log10(mG10h_y(mG10h_y>=0.1)),20,'ko','filled')
% plot(G10h_x,-log10(G10h_y),'bo')
set(gca,'LineWidth',1.25,'FontName','Arial','FontSize',12);
xlabel('log_2(Fold change)'); ylabel('-log_1_0(FDR)');
text(-7.9,62,'d','FontName','Arial','FontSize',14)
axis([-6 6 0 60])
xticks([-6:3:6]);
yticklabels = get(gca, 'YTickLabel');
yticklabels{2}='';yticklabels{4}='';yticklabels{6}='';
set(gca, 'YTickLabel', yticklabels);
% lg=text(3.1,56,'APAP (10h)');
lg.FontSize = 11;lg.FontName = 'arial';


