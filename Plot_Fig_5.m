changeCobraSolver('glpk'); % or glpk
% % Load COBRA models
% remove "insertpath" and  put current folder location of the files
% downloaded for manuscript
addpath '"Insertpath"'
load iRno_v2.mat;
model=rno_cobra;
% Find total number of unique metabolites in the model
UkegID = unique(model.metKEGGID,'sorted');
MkegID = UkegID(2:1701);
%% matching metabolites from plasma (blood) sample to the model
addpath '"Insertpath"'
% Kegg ids in the data
[~,T,Kraw] = xlsread('TableS3','Pathway Heat Map','H9:H577');
[v5h,~,Kraw5] = xlsread('TableS3','Pathway Heat Map','Q9:Q577');
[v10h,~,Kraw10] = xlsread('TableS3','Pathway Heat Map','R9:R577');
[p5h,~,Krawp5] = xlsread('TableS3','Pathway Heat Map','AD9:AD577');
[p10h,~,Krawp10] = xlsread('TableS3','Pathway Heat Map','AF9:AF577');
[q5h,~,Krawq5] = xlsread('TableS3','Pathway Heat Map','AE9:AE577');
[q10h,~,Krawq10] = xlsread('TableS3','Pathway Heat Map','AG9:AG577');
[~,N,Nraw] = xlsread('TableS3','Pathway Heat Map','E9:E577');

% volcano plot showing metabolite changes for 5 h and 10h
figure (5); set(figure(5),'Units','inches','Position',[0.5 0.5 9.8 7.0]);
subplot(2,2,1); subplot('position',[0.08 0.58 0.4 0.38]);
hold on; box off; grid off;
ve5h = log2(v5h); %y5h = -log10(p5h);
scatter(ve5h(q5h<0.1 & ve5h>0),-log10(q5h(q5h<0.1 & ve5h>0)),20,'ro','filled');hold on;
scatter(ve5h(q5h<0.1 & ve5h<0),-log10(q5h(q5h<0.1 & ve5h<0)),20,'go','filled');hold on;
scatter(ve5h(q5h>=0.1),-log10(q5h (q5h>=0.1)),20,'ko','filled');
set(gca,'LineWidth',1.5,'FontSize',12,'FontName','arial');
xlabel('log_2(Fold change)','FontName','arial'); ylabel('-log_1_0(FDR)')
text(-7.6,12.6,'a','FontName','Arial','FontSize',14)
xlim([-6,6]); xticks(-6:3:6);ylim([0,12]);yticks(0:3:12)
lg=text(3.35,12.6,'APAP (5 h)');
lg.FontSize = 11;lg.FontName = 'arial';

subplot(2,2,2); subplot('position',[0.57 0.58 0.4 0.38]);
hold on; box off; grid off;
ve10h = log2(v10h); 
scatter(ve10h(q10h<0.1 & ve10h>0),-log10(q10h(q10h<0.1 & ve10h>0)),20,'ro','filled');hold on;
scatter(ve10h(q10h<0.1 & ve10h<0),-log10(q10h(q10h<0.1 & ve10h<0)),20,'go','filled');hold on;
scatter(ve10h(q10h>=0.1),-log10(q10h(q10h>=0.1)),20,'ko','filled');
set(gca,'LineWidth',1.5,'FontSize',12,'FontName','arial');
xlabel('log_2(Fold change)','FontName','arial'); ylabel('-log_1_0(FDR)')
text(-7.6,12.6,'b','FontName','Arial','FontSize',14)
xlim([-6,6]); xticks(-6:3:6);ylim([0,12]);yticks(0:3:12)
lg=text(3.35,12.6,'APAP (10 h)');
lg.FontSize = 11;lg.FontName = 'arial';

%% Corresponding Biochemical name
[~,Mname] = xlsread('TableS3','Pathway Heat Map','E9:E577');
% find unique metabolite ids removing NaN's in the Kraw
UKraw = unique(T,'sorted');
DKraw = UKraw(2:end);
% compare matching kegg id numbers in the model and data
[i,j]=ismember(DKraw,MkegID);
% number of metabolites in the data that are in the model
n_m = find(i~=0);
% Find kegg ids for matched metabolites between data and model
M_mn = j(n_m);
M_m = MkegID(M_mn);
% find the metaboites names for the above keggs ids in the model
[Met_keg] = [model.metNames model.metKEGGID];
[p,q] = ismember(M_m, model.metKEGGID);
% Names of matched metabolites in the model
met_names_model = model.metNames(q);
% Metnames from the data
[p1,q1]=ismember(M_m,T);
met_names_data=Mname(q1);
% Load the metabolites that were mapped to the iRno model from TableS3 based on KEGG ID
% annotation
All_Data10h=readtable('"Insertpath"\Metabolites_mapped_10h_data.txt','Delimiter','\t');
All_Data5h=readtable('"Insertpath"\Metabolites_mapped_5h_data.txt','Delimiter','\t');
M10 = All_Data10h{:,1};
M5 = All_Data5h{:,1};
[p4,q4]=ismember(M5,M10);
M5nt10 = find(p4~=1);
% Model exchangeable metabolites that are in the data for both 5 h and 10h
all_mets_exch = [M10;M5(M5nt10)];
% find metabolites that are not in the list that mapped to data by keggID comparison with metabolites
% that were identified based on names from different databases in % all_mets_exch
[p2,q2]=ismember(all_mets_exch,met_names_model);
Nim = find(p2~=1);
Mtn = all_mets_exch(Nim);
Mtn(1)={'S-1-pyrroline-5-carboxylate'};% kegg id exists in data but not in the model
Mtn(5) = {'2-hydroxybutyrate/2-hydroxyisobutyrate'};
Mtn(7) = {'taurohyodeoxycholic acid'};
Mtn(8)={'arabitol/xylitol'};
Mtn(9)={'oleate/vaccenate (18:1)'};
Mtn(10) = {'oleoylcarnitine'};
Mtn(11) = {'tauro-alpha-muricholate'};
Mtn(12) = {'tauro-beta-muricholate'};
Mtn(13) = {'glycohyodeoxycholate'};
Mtn(15) = {'margarate (17:0)'};
Mtn(16) = {'10-heptadecenoate (17:1n7)'};
Mtn(18) = {'mead acid (20:3n9)'};
Mtn(19) = {'arabonate/xylonate'};
Mtn(20) = {'3-methyl-2-oxovalerate'};
Mtn(21)={'alpha-hydroxyisocaproate'};
Mtn(22) = {'3-phosphoglycerate'};
Mtn(24) = {'arabitol/xylitol'};
% names matching between model metabolites and data
[p3,q3]=ismember(Mtn,Mname);
q1 = [q1;q3(q3~=0)];
Kid = T(q1);v5h_m = v5h(q1); p5h_m = p5h(q1);v10h_m = v10h(q1);p10h_m = p10h(q1);
q5h_m = q5h(q1);q10h_m = q10h(q1);
met_names_data=Mname(q1);
% Metabolites that are in the model and the corresponding fold change
% values and p values based on KEGG ID comparision
matched_mets = [met_names_data(:,1) Kid(:,1) num2cell(v5h_m(:,1)) num2cell(q5h_m(:,1)) num2cell(v10h_m(:,1)) num2cell(q10h_m(:,1))];
% figure(2); hold on;% volcano plot showing metabolite changes for 5 h and 10h that are mapped to the model
subplot(2,2,3); subplot('position',[0.08 0.1 0.4 0.38]);
hold on; box off; grid off;
v5hx_m = log2(v5h_m);
scatter(v5hx_m(q5h_m<0.1 & v5hx_m>0),-log10(q5h_m(q5h_m<0.1 & v5hx_m>0)),20,'ro','filled');hold on;
scatter(v5hx_m(q5h_m<0.1 & v5hx_m<0),-log10(q5h_m(q5h_m<0.1 & v5hx_m<0)),20,'go','filled');hold on;
scatter(v5hx_m(q5h_m>=0.1),-log10(q5h_m(q5h_m>=0.1)),20,'ko','filled');hold on;
set(gca,'LineWidth',1.5,'FontSize',12,'FontName','arial');
xlabel('log_2(Fold change)','FontName','arial'); ylabel('-log_1_0(FDR)')
text(-7.6,9.6,'c','FontName','Arial','FontSize',14)
xlim([-6,6]); ylim([0 9]); xticks(-6:3:6);yticks(0:3:9)

subplot(2,2,4); subplot('position',[0.57 0.1 0.4 0.38]);
hold on; box off; grid off;
v10hx_m = log2(v10h_m);
scatter(v10hx_m(q10h_m<0.1 & v10hx_m>0),-log10(q10h_m(q10h_m<0.1 & v10hx_m>0)),20,'ro','filled');hold on;
scatter(v10hx_m(q10h_m<0.1 & v10hx_m<0),-log10(q10h_m(q10h_m<0.1 & v10hx_m<0)),20,'go','filled');hold on;
scatter(v10hx_m(q10h_m>=0.1),-log10(q10h_m(q10h_m>=0.1)),20,'ko','filled');hold on;
set(gca,'LineWidth',1.5,'FontSize',12,'FontName','arial');
xlabel('log_2(Fold change)','FontName','arial'); ylabel('-log_1_0(FDR)')
text(-7.6,9.6,'d','FontName','Arial','FontSize',14)
xlim([-6,6]); ylim([0 9]); xticks(-6:3:6);yticks(0:3:9)


