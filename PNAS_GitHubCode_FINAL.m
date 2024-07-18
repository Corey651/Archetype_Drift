%% C. Elegans Analysis

%% N-NMF Analysis
T=readtable('Petrascheck_Cleaned_Data.xlsx');

% Analysis based on sample means
Data=T{:,15:18};
Names=T.GeneID;
M=mean(Data,2);
Data2=Data;
Data=Data./M;
V=var(Data,0,2);
Ind=find(V>=.1&V<=.2);

Gene_Labels=Names(Ind);
DD=Data;


DD=DD(Ind,:);

% Figure S1

[u,s,v]=svd(DD);
figure
plot(diag(s))
xlabel('Number of Components')
ylabel('Eigenvalues')

DD=DD./sum(DD,1);

[W,H]=NMF(DD,2);

% Plot the Archetypes (Figure 1)
Time=[1,3,5,10];
difftime = Time - 1;

figure
axes1 = axes('Position',...
    [0.127395833333333 0.178141592920354 0.775 0.774513274336283]);
hold(axes1,'on');
bar1 = bar(Time,H,'BarLayout','stacked');
set(bar1(2),...
    'FaceColor',[0.466666666666667 0.674509803921569 0.188235294117647]);
ylabel('Archetype Score');
xlabel('Days');
box(axes1,'on');
hold(axes1,'off');


% Calculate Correlations of Genes with Archetype 2
CC=corrcoef([Data2;H]');
CC=CC(1:20537,20539);

% Plot Correlations as a Histogram for All Genes (Figure 2A)
figure
histogram(CC,1000)
xlabel('Correlation')
ylabel('Count')

Corr_Table=table('Size',[20537 2],'VariableTypes',{'string','double'},'VariableNames',{'Gene_Name','Correlation'});
Corr_Table.Gene_Name=Names;
Corr_Table.Correlation=CC;
%writetable(Corr_Table,'Correlation_Table.csv')
hist(Corr_Table.Correlation, 40)
hold on
xline(0.9)

%% Decreasing Genes 
% Plot Decreasing Archetype Centers (Figure 2B)
figure 
plot(Time,Data(Corr_Table.Correlation<=-.90,:)')

% Plot Mean Percent Change in Decreasing Genes
figure
decreasinggenes = Data(Corr_Table.Correlation<=-.90,:);
meandecreasing = mean(decreasinggenes);
percentchangedec = 100*(meandecreasing - meandecreasing(:, 1))/meandecreasing(:, 1);
scatter(difftime, percentchangedec)

% Fit Model to Percent Change in Decreasing Genes
start2 = [1, 1];
lower2 = [0,0];
ptime = [0:0.1:9];
fitfun = fittype( @(g,K,x)  (g*K)./(K+ 0.1963*x) - g);
[fitted_curve,gof] = fit(difftime.',percentchangedec.',fitfun,'StartPoint',start2, 'Lower', lower2);

% Plot results
hold on
plot(ptime,fitted_curve(ptime),'LineWidth',1.5);
scatter(difftime, percentchangedec, 30)
% Calculates R^2 and Stores in Vector
rsquares = gof.rsquare;
% Stores “g” and “K” values in Vector (with “g” in first column)
coeffvals = [coeffvalues(fitted_curve)];


% % Individual Fits Decreasing Genes (Figure S2B)
% % This code takes a while to run.
% ptime = [0:0.1:9];
% t2 = [0, 2, 4, 9]';
% individualpercentchangedecreasing = 100*(decreasinggenes - decreasinggenes(:, 1))./decreasinggenes(:, 1);
% individualchangedecreasing = (decreasinggenes - decreasinggenes(:, 1));
% colors = [0 0.447 0.741; 0.929 0.694 0.125; 0.466 0.674 0.188; 0.635 0.078 0.184; 0.85 0.325 0.098; 0.494 0.184 0.556; 0.301 0.745 0.933];
% % Define Start points, fit-function
% colorindex = 1;
% start2 = [1, 1];
% lower2 = [0,0];
% fitfun = fittype( @(g,K,x)  (g*K)./(K+ 0.1963*x) - g);
%
% figure
% for c = 1:1859
% %Fit curve
% dummy = (individualchangedecreasing(c,:))';
% [fitted_curve,gof] = fit(t2,dummy,fitfun,'StartPoint',start2, 'Lower', lower2);
% % Plot results
% hold on
% plot(ptime,fitted_curve(ptime),'LineWidth',1.5);
% scatter(t2, dummy, 30, colors(colorindex, :),  'filled')
% % Calculates R^2 and Stores in Vector
% rsquares(c, 1) = gof.rsquare;
% % Stores “g” and “K” values in Vector (with “g” in first column)
% coeffvals(c, 1:2) = [coeffvalues(fitted_curve)];
% colorindex = colorindex + 1;
% if colorindex == 8
% colorindex = 1;
% end
% end
% % Histogram of Rsquares (Figure S2B)
% figure
% edges = linspace(0, 1, 21); % Create 20 bins.
% % Plot the histogram.
% histogram(rsquares, 'BinEdges',edges);
% % Fancy up the graph.
% xlim([0, 1]);
% xlabel('R^2', 'FontSize', 14);
% ylabel('Bin Count', 'FontSize', 14);


%% Increasing Genes
% Plot Increasing Archetype Centers (Figure 2B)
figure
plot(Time,Data(Corr_Table.Correlation>=.90,:)') 

% Increasing Genes Mean Percent Change (Figure 3A)
figure
increasinggenes = Data(Corr_Table.Correlation>=.90,:);
meanincreasing = mean(increasinggenes);
percentchangeinc = 100*(meanincreasing - meanincreasing(:, 1))/meanincreasing(:, 1);
scatter(difftime, percentchangeinc)

% Fit Equation to Mean Percent Change of Increasing Genes
start2 = [1];
lower2 = [0];
ptime = [0:0.1:9];
fitfun = fittype( @(a,x)  (a*x));
[fitted_curve,gof] = fit(difftime.',percentchangeinc.',fitfun,'StartPoint',start2, 'Lower', lower2);
% Plot results
hold on
plot(ptime,fitted_curve(ptime),'LineWidth',1.5);
scatter(difftime, percentchangeinc, 30)
% Calculates R^2 and Stores in Vector
rsquares = gof.rsquare;
% Stores coeff values in Vector
coeffvals = [coeffvalues(fitted_curve)];
hold off


% % Individual Fits Increasing Genes (Figure S2A)
% % This code takes a while to run. 
% ptime = [0:0.1:9];
% t2 = [0, 2, 4, 9]';
% individualpercentchangeincreasing = rmmissing(100*(increasinggenes - increasinggenes(:, 1))./increasinggenes(:, 1));
% individualchangeincreasing = rmmissing((increasinggenes - increasinggenes(:, 1))); % Removes ~800 genes which have a NaN divide by 0 error
% colors = [0 0.447 0.741; 0.929 0.694 0.125; 0.466 0.674 0.188; 0.635 0.078 0.184; 0.85 0.325 0.098; 0.494 0.184 0.556; 0.301 0.745 0.933];
% colorindex = 1
% start1 = 1;
% lower1 = 0;
% fitfun = fittype( @(b, x) b*x);
% figure
% for c = 1:3006
% % Define Start points, fit-function and fit curve
% dummy = (individualchangeincreasing(c,:))';
% [fitted_curve,gof] = fit(t2,dummy,fitfun,'StartPoint',start1,'Lower', lower1);
% % Plot results
% hold on
% plot(ptime,fitted_curve(ptime),'LineWidth',1.5);
% scatter(t2, dummy, 30, colors(colorindex, :),  'filled')
% colorindex = colorindex + 1;
% % Calculates R^2 and Stores in Vector
% rsquares(c, 1) = gof.rsquare;
% % Stores “b” value in Vector 
% coeffvals(c, 1) =[coeffvalues(fitted_curve)];
% if colorindex == 8
% 	colorindex = 1;
% end
% end
% % Histogram of Rsquares (Figure S2A)
% figure
% edges = linspace(0, 1, 21); % Create 20 bins.
% % Plot the histogram.
% histogram(rsquares, 'BinEdges',edges);
% % Fancy up the graph.
% xlim([0, 1]);
% xlabel('R^2', 'FontSize', 14);
% ylabel('Bin Count', 'FontSize', 14);


%% Human Fibroblast Work
%% Archetype Analysis 
TPic=readtable('Picard_Human_Fibroblasts.csv');
Genes_Pic=TPic.Var1;
Subj=TPic.Properties.VariableNames;
Subj=Subj(2:end);
Subj=extractBetween(Subj,'x','_');

TPic=TPic{:,2:end};

T_Data=TPic;

T_Meta=readtable('Picard_Human_Fibroblasts_Metadata.csv');
RNA_Order=T_Meta.RNAseq_sampleID;
[~,ia,ib]=intersect(RNA_Order,str2double(Subj));

RNA_Order=RNA_Order(ia);
T_Meta=T_Meta(ia,:);
T_Data=T_Data(:,ib);


Grps=unique(T_Meta.Cell_line_group);

Grp_Cell=cell(1,length(Grps));
for i=1:length(Grp_Cell)
    Grp_Cell{i}=find(strcmp(T_Meta.Cell_line_group,Grps(i)));
end

Grp_size=zeros(1,length(Grps));
for i=1:length(Grps)
    Grp_size(i)=length(Grp_Cell{i});
end


% Train archetypes on everything
Data_Control_1=T_Data(:,[Grp_Cell{1}]);
T_Meta_Control_1=T_Meta([Grp_Cell{1}],:);
Time_1=T_Meta_Control_1.Days_grown_Udays;

Data_Control_2=T_Data(:,[Grp_Cell{15}]);
T_Meta_Control_2=T_Meta([Grp_Cell{15}],:);
Time_2=T_Meta_Control_2.Days_grown_Udays;

Data_Control_3=T_Data(:,[Grp_Cell{29}]);
T_Meta_Control_3=T_Meta([Grp_Cell{29}],:);
Time_3=T_Meta_Control_3.Days_grown_Udays;



Data_Control=(Data_Control_1+Data_Control_2)/2;
Time=(Time_1+Time_2)/2;

M=mean(Data_Control,2);
Data_Control=Data_Control./M; 

V=var(Data_Control,0,2);
Ind=find(V>=.4&V<=.6);



Data_Surf1_1=T_Data(:,[Grp_Cell{40}]);
T_Meta_Surf1_1=T_Meta([Grp_Cell{40}],:);
Time_M1=T_Meta_Surf1_1.Days_grown_Udays;

Data_Surf1_2=T_Data(:,[Grp_Cell{42}]);
T_Meta_Surf1_2=T_Meta([Grp_Cell{42}],:);
Time_M2=T_Meta_Surf1_2.Days_grown_Udays;


Data_Control_1=Data_Control_1./M;
Data_Control_2=Data_Control_2./M;
Data_Control_3=Data_Control_3./M;

Data_Surf1_1=Data_Surf1_1./M; 
Data_Surf1_2=Data_Surf1_2./M; 


Data_Control_1=Data_Control_1(Ind,:);
Data_Control_2=Data_Control_2(Ind,:);
Data_Control_3=Data_Control_3(Ind,:);

Data_Surf1_1=Data_Surf1_1(Ind,:);
Data_Surf1_2=Data_Surf1_2(Ind,:);


DC1=Data_Control_1./sum(Data_Control_1,1);
DC2=Data_Control_2./sum(Data_Control_2,1);
DC3=Data_Control_3./sum(Data_Control_3,1);

DM1=Data_Surf1_1./sum(Data_Surf1_1,1);
DM2=Data_Surf1_2./sum(Data_Surf1_2,1);

[W,H]=NMF((DC1+DC2)/2,2);

[u,s,v]=svd((DC1+DC2)/2);  % plot(diag(s)) to see that 1 component captures most of the variance

HC3=NMF_New_Weights(DC3,W,2);

HM=NMF_New_Weights([DM1+DM2]/2,W,2);  % Use this to calculate archetype scores from a trained model

% Archetype Plot (Figure 4)
figure
yyaxis left
plot(Time,H(1,:))
hold on
plot(Time_3,HC3(1,:))
plot((Time_M1+Time_M2)/2,HM(1,:))
ylabel('Archetype 1')

yyaxis right
plot(Time,H(2,:))
hold on
plot(Time_3,HC3(2,:))
plot((Time_M1+Time_M2)/2,HM(2,:))
ylabel('Archetype 2')

set(gcf, 'color', 'w')
xlabel('Days Grown')
legend
hold off

%% Retrain Archetype on just Control Datasets
T=readtable('Picard_Human_Fibroblasts_Cleaned_Controls');

Data=T{:,2:12};
Names=T.Gene;
M=mean(Data,2);
Data2=Data;
Data=Data./M;
V=var(Data,0,2);
Ind=find(V>=.1&V<=.2);
Gene_Labels=Names(Ind);
DD=Data;

DD=DD(Ind,:);
DD=DD./sum(DD,1);
[W,H]=NMF(DD,2);

Time = [21.21, 37.33, 52.96, 68.08, 83.08, 98.08, 115, 134.96, 157.92, 180.88, 202];
difftime = Time - 21.21;

CC=corrcoef([Data2;H]');
CC=CC(1:26294,26296);

% Plot Correlation Histogram (Figure 5A)
figure
histogram(CC,1000)
xlabel('Correlation')
ylabel('Count')

Corr_Table=table('Size',[26294 2],'VariableTypes',{'string','double'},'VariableNames',{'Gene_Name','Correlation'});
Corr_Table.Gene_Name=Names;
Corr_Table.Correlation=CC;
% writetable(Corr_Table,'Correlation_Table_HC.csv')

%% Increasing Genes
% Plot Increasing Archetype Centers (Figure 5B)
figure
plot(Time,Data(Corr_Table.Correlation>=.90,:)')

% Figure S3A
% Plot Percent Change in Increasing Archetype
figure
increasinggenes = Data(Corr_Table.Correlation>=.90,:);
meanincreasing = mean(increasinggenes);
percentchangeinc = 100*(meanincreasing - meanincreasing(:, 1))/meanincreasing(:, 1);
scatter(difftime, percentchangeinc)
set(gcf, 'color', 'w')
xlabel('Days')
ylabel('Percent Change in Mean Gene Expression (CPM)')

% Hill Fit to Increasing Genes
% Define Start points, fit-function
start2 = [1, 1, 1];
lower2 = [0,0, 0];
ptime = [0:0.5:180.79];
fitfun = fittype( @(a,b,n,x)  (b*x.^n)./(a^n + x.^n));

%Fit curve
[fitted_curve,gof] = fit(difftime.',percentchangeinc.',fitfun,'StartPoint',start2, 'Lower', lower2);
% Plot results
hold on
plot(ptime,fitted_curve(ptime),'LineWidth',1.5);
scatter(difftime, percentchangeinc, 30)
% Calculates R^2 and Stores in Vector
rsquares = gof.rsquare;


%% Decreasing Genes
% Plot Decreasing Archetype Centers (Figure 5B)
figure
plot(Time,Data(Corr_Table.Correlation<=-.90,:)')

% Figure S3B
% Plot Percent Change in Decreasing Archetype
figure
decreasinggenes = Data(Corr_Table.Correlation<=-.90,:);
meandecreasing = mean(decreasinggenes);
percentchangedec = 100*(meandecreasing - meandecreasing(:, 1))/meandecreasing(:, 1);
scatter(difftime, percentchangedec)
set(gcf, 'color', 'w')
xlabel('Days')
ylabel('Percent Change in Mean Gene Expression (CPM)')
% Decreasing Genes Fit
start2 = [1, 1];
lower2 = [0,0];
ptime = [0:0.5:180.79];
fitfun = fittype( @(g,K,x)  (g*K)./(K+ 0.1963*x) - g);

[fitted_curve,gof] = fit(difftime.',percentchangedec.',fitfun,'StartPoint',start2, 'Lower', lower2);
% Plot results
hold on
plot(ptime,fitted_curve(ptime),'LineWidth',1.5);
scatter(difftime, percentchangedec, 30)
% Calculates R^2 and Stores in Vector
rsquares = gof.rsquare;
% Stores “g” and “K” values in Vector (with “g” in first column)
coeffvals = [coeffvalues(fitted_curve)];



