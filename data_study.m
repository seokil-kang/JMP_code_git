% data_index: vector of integers with the following order
% 1 2 3 4 5 6 7 8 9 10 11 12 13 14
% y c i w n b s g z p  R  Rb x s(accounting)
% t1: data starting quarterly date e.g. 'yyyy-Q' or '1960-1'
% t2: data ending quarterly date e.g. 'yyyy-Q' or '2006-4''
% quarterly date means the last digit can take value 1~4
% differenced: 1: data to be stationary 2: data to be level

% house keeping
clear; close all; clc;

% change directory to data folder
cd data

% specify full sample period
t1 = '1960-1';
t2 = '2006-4';

% taking difference to nonstationary observable?
differenced = 1;

% which data to select?
data_index = linspace(1,14,14)';

% benchmark data selection
%data_index = [1 3 4 5 6 7 10 11 12];

% HPC does not have the X-13 filter add-on
isHPC = 0;

% quarterly data range from 1954Q3 to 2021Q1
data_quarterly = readmatrix('fredgraphQ.xls','range','B49:I315');
date_quarterly_full = [datetime(1954,7,1,'format','yyyy-Q'):calmonths(3):datetime(2021,1,1,'format','yyyy-Q')]';
% Gross Domestic Product, Billions of Dollars, Quarterly, Seasonally Adjusted Annual Rate
GDP = data_quarterly(:,1);
% Personal Consumption Expenditures, Billions of Dollars, Quarterly, Seasonally Adjusted Annual Rate
PCE = data_quarterly(:,2);
% Fixed Private Investment, Billions of Dollars, Quarterly, Seasonally Adjusted Annual Rate
FPI = data_quarterly(:,3);
% Gross Domestic Product: Implicit Price Deflator, Index 2012=100, Quarterly, Seasonally Adjusted
GDPDEF = data_quarterly(:,4);
% Nonfarm Business Sector: Compensation Per Hour, Index 2012=100, Quarterly, Seasonally Adjusted
COMPNFB = data_quarterly(:,5);
% Federal Consumption Expenditures and Gross Investment, Billions of Dollars, Quarterly, Seasonally Adjusted Annual Rate
GCI = data_quarterly(:,6);
% Federal government current transfer payments: Government social benefits: to persons, Billions of Dollars, Quarterly, Seasonally Adjusted Annual Rate
GTRSF = data_quarterly(:,7);
% Nonfarm Business Sector: Average Weekly Hours, Index 2012=100, Quarterly, Seasonally Adjusted
AWH = data_quarterly(:,8);
% Nominal Potential Gross Domestic Product, Billions of Dollars, Quarterly, Not Seasonally Adjusted
PGDP = readmatrix('NGDPPOT.xls','range','B34:B300');
% Unemployment Rate, Percent, Quarterly, Seasonally Adjusted
UMP = readmatrix('UNRATE.xls','range','B39:B304');

% monthly data range from 1954M7 to 2021M3
data_monthly = readmatrix('fredgraphM.xls','range','B168:H968');
date_monthly = [datetime(1954,7,1,'format','yyyy-M'):calmonths(1):datetime(2021,3,1,'format','yyyy-M')]';

% Effective Federal Funds Rate, Percent, Monthly, Not Seasonally Adjusted
FFR_m = data_monthly(:,1);
% Employment Level, Thousands of Persons, Monthly, Seasonally Adjusted
EMP_m = data_monthly(:,2);
% Normalize it at 2012
EMP_m = EMP_m / mean(EMP_m(find(date_monthly == '2012-1'):find(date_monthly == '2012-12')));
% Population Level, Thousands of Persons, Monthly, Not Seasonally Adjusted
POP_m = data_monthly(:,3);
% Normalize it at 2012
POP_m = POP_m / mean(POP_m(find(date_monthly == '2012-1'):find(date_monthly == '2012-12')));
% Market Value of Privately Held Gross Federal Debt, Billions of Dollars, Monthly, Not Seasonally Adjusted
MKVDDF_m = data_monthly(:,4);
% 5-Year Treasury Constant Maturity Rate, Percent, Monthly, Not Seasonally Adjusted
CMR5_m = data_monthly(:,5);
% 10-Year Treasury Constant Maturity Rate, Percent, Monthly, Not Seasonally Adjusted
CMR10_m = data_monthly(:,6);
% 3-Year Treasury Constant Maturity Rate, Percent, Monthly, Not Seasonally Adjusted
CMR3_m = data_monthly(:,7);

% Effective Federal Funds Rate, Percent, Monthly, Not Seasonally Adjusted
FFR = mean(reshape(FFR_m,3,length(date_quarterly_full)))'/4;
% Employment Level, Thousands of Persons, Monthly, Seasonally Adjusted
EMP = mean(reshape(EMP_m,3,length(date_quarterly_full)))';
% Population Level, Thousands of Persons, Monthly, Not Seasonally Adjusted
POP = mean(reshape(POP_m,3,length(date_quarterly_full)))';
% Market Value of Privately Held Gross Federal Debt, Billions of Dollars, Monthly, Not Seasonally Adjusted
MKVDDF = mean(reshape(MKVDDF_m,3,length(date_quarterly_full)))';
% 3-Year Treasury Constant Maturity Rate, Percent, Monthly, Not Seasonally Adjusted
CMR3 = mean(reshape(CMR3_m,3,length(date_quarterly_full)))'/4;
% 5-Year Treasury Constant Maturity Rate, Percent, Monthly, Not Seasonally Adjusted
CMR5 = mean(reshape(CMR5_m,3,length(date_quarterly_full)))'/4;
% 10-Year Treasury Constant Maturity Rate, Percent, Monthly, Not Seasonally Adjusted
CMR10 = mean(reshape(CMR10_m,3,length(date_quarterly_full)))'/4;

% US treasury holding period return by Hall, Payne and Sargent(2018), Percent, Monthly, Not Seasonally Adjusted
HPR_m = readmatrix('US_holding_period_return_1790_2021.xlsx','range','C1975:C2775')*1e2;
% US treasury holding period return US treasury debt data by Hall, Payne and Sargent(2018) Billons of Dollars, Monthly, Not Seasonally Adjusted
MKVD_m = readmatrix('US_Treasury_Debt_1776-2021.xlsx','range','I2135:I2935')/1e9;

% quarterlized HPR and MKVD
HPR = prod(reshape(1+HPR_m/1e2,3,length(date_quarterly_full)))';
HPR = (HPR-1)*1e2;
MKVD3 = (reshape(MKVD_m,3,length(date_quarterly_full)))';
%MKVD = MKVD3(:,3); % Cochrane uses this way
MKVD = mean(MKVD3,2);
% debt to GDP ratio
BYratio = MKVD ./ GDP;
% growth rate with consumption
growth_rate = diff(1e2*log(GDP./GDPDEF));
growth_rate_per_capita = diff(1e2*log(GDP./GDPDEF./POP));

% primary surplus accounting method
% governemt part
date_ps = [datetime(1959,7,1,'format','yyyy-Q'):calmonths(3):datetime(2018,10,1,'format','yyyy-Q')]';
% total receipts: Table 3.2 line 40
rcpt = [readmatrix('NIPA_table_3_2.xls','sheet','Sheet0','range','BA48:IV48')' ; readmatrix('NIPA_table_3_2.xls','sheet','Sheet1','range','C48:AJ48')'];

% total expenditure: Table 3.2 line 43
expd = [readmatrix('NIPA_table_3_2.xls','sheet','Sheet0','range','BA51:IV51')' ; readmatrix('NIPA_table_3_2.xls','sheet','Sheet1','range','C51:AJ51')'];

% current expenditure: Table 3.2 line 44
expd_crnt = [readmatrix('NIPA_table_3_2.xls','sheet','Sheet0','range','BA52:IV52')' ; readmatrix('NIPA_table_3_2.xls','sheet','Sheet1','range','C52:AJ52')'];

% gross government investment: Table 3.2 line 45
govt_inv = [readmatrix('NIPA_table_3_2.xls','sheet','Sheet0','range','BA53:IV53')' ; readmatrix('NIPA_table_3_2.xls','sheet','Sheet1','range','C53:AJ53')'];

% capital transfer payments: Table 3.2 line 46
k_trans_pay = [readmatrix('NIPA_table_3_2.xls','sheet','Sheet0','range','BA54:IV54')' ; readmatrix('NIPA_table_3_2.xls','sheet','Sheet1','range','C54:AJ54')'];

% (Less:) Consumption of fixed capital: Table 3.2 line 48
cons_fix_k = [readmatrix('NIPA_table_3_2.xls','sheet','Sheet0','range','BA56:IV56')' ; readmatrix('NIPA_table_3_2.xls','sheet','Sheet1','range','C56:AJ56')'];

% interest payments: Table 3.2 line 33
int_pay = [readmatrix('NIPA_table_3_2.xls','sheet','Sheet0','range','BA40:IV40')' ; readmatrix('NIPA_table_3_2.xls','sheet','Sheet1','range','C40:AJ40')'];

% interest receipts: Table 3.2 line 14
int_rcpt = [readmatrix('NIPA_table_3_2.xls','sheet','Sheet0','range','BA21:IV21')' ; readmatrix('NIPA_table_3_2.xls','sheet','Sheet1','range','C21:AJ21')'];

% naive interpolation
int_rcpt(2) = int_rcpt(3);
int_rcpt(1) = int_rcpt(2);

% Fed employee pension and insurance fund trans => interest accrued: Table 3.18A and 3.18B line 22
pension_pay = [readmatrix('NIPA_table_3_18A.xls','sheet','Sheet0','range','C31:AJ31')' ; readmatrix('NIPA_table_3_18B.xls','sheet','Sheet0','range','C31:GX31')'];

% compute primary surplus (Hall&Sargent, 2011)
expd_adjst = expd - (int_pay - int_rcpt + pension_pay);
lnps_account = (rcpt - expd_adjst) ./ GDPDEF(find(date_quarterly_full == date_ps(1)):find(date_quarterly_full == date_ps(end)));

% Primary surplus to GDP ratio
SYratio_identity = -diff(1e2*log(BYratio)) + HPR(2:end) - diff(1e2*log(GDPDEF)) - growth_rate;

% calibrated time discount rate
betat = .99;

% log of primary surplus from gov't identity
lnps_identity_NSA = -diff(1e2*log(MKVD./ GDPDEF ./ POP)) + 1e2*log(1+HPR(2:end)/1e2) - diff(1e2*log(GDPDEF));

% compare the data with Cochrane RED 2021
load surplus_Cochrane
surplus_Cochrane = [surplus_identity_q(31:end); lnps_identity_NSA(259:end)];

% plot the primary surplus series
if 0
    figure
    hold on
    p1 = plot(date_quarterly_full(2:end),SYratio_identity,'linewidth',1.5);
    p2 = plot(date_quarterly_full(2:end),lnps_identity_NSA,'linewidth',1.5);
    p3 = plot(datec,surplus_identity_q*1e2,'g','linewidth',1);
    hold off
    recessionplot
    legend([p1 p2 p3],'SYratio','lnPS', 'Cochrane','location','southoutside','orientation','horizontal')
    legend boxoff
end

% skip the X-13 filter for HPC case
if isHPC
    load seasonally_adjusted_primary_surplus
else
    % primary surplus suffers seasonality a lot...
    % apply X-13 filter
    spec = makespec('FLOW','TRAMO','SEATS','series','comptype','add','CONSTANT');
    X13_ps = x13([datenum(date_quarterly_full(2:end)),lnps_identity_NSA], x13spec(spec,'series','name','primary_surplus'));
    lnps_identity = X13_ps.s11.s11;
    save('seasonally_adjusted_primary_surplus.mat','lnps_identity')
end

% go back to the original directory
cd ..

% Construct log real per-capita values
lny = 1e2 * log(GDP ./ GDPDEF ./ POP);
lnc = 1e2 * log(PCE ./ GDPDEF ./ POP);
lni = 1e2 * log(FPI ./ GDPDEF ./ POP);
lnw = 1e2 * log(COMPNFB ./ GDPDEF);
lnn = 1e2 * log(AWH .* EMP ./ POP);
lnp = 1e2 * log(GDPDEF);
lnb = 1e2 * log(MKVD ./ GDPDEF ./ POP);
lng = 1e2 * log(GCI ./ GDPDEF ./ POP);
lnz = 1e2 * log(GTRSF ./ GDPDEF ./ POP);
lnyf = 1e2 * log(PGDP ./ GDPDEF ./ POP);

% define output gap
lnx = lny-lnyf;

% take out the first point due to differencing
date_quarterly = date_quarterly_full(2:end);

% primary surplus accounting misses values in pre 1959 and post 2018
lnps_account = [nan((find(date_quarterly == date_ps(1))-1),1); lnps_account; nan(length(date_quarterly)-find(date_quarterly == date_ps(end)),1)];

% full sample in level
y_full = [lny(2:end) lnc(2:end) lni(2:end) lnw(2:end) lnn(2:end) lnb(2:end) lnps_identity lng(2:end) lnz(2:end) diff(lnp) FFR(2:end) HPR(2:end) lnx(2:end) lnps_account];
% full sample differenced
dy_full = [diff(lny) diff(lnc) diff(lni) diff(lnw) lnn(2:end) diff(lnb) lnps_identity diff(lng) diff(lnz) diff(lnp) FFR(2:end) HPR(2:end) lnx(2:end) lnps_account];
% record differenced variables
differenced_index_full = [1 2 3 4 6 8 9];
% full sample data label
y_nm_full = ["output";"consumption";"investment";"wage";"labor";"debt";"surplus";"govt spending";"transfer";"inflation";"interest rate";"holding return";"output gap";"surplus"];

% plot the full sample
if 0
    for j = 1:length(y_nm_full)
        subplot(7,2,j)
        plot(date_quarterly,dy_full(:,j),'.-','linewidth',1.5)
        recessionplot
        axis tight
        title(y_nm_full(j))
    end
end

% differenced_set
differenced_index = [];

% choose data differencing
if differenced == 0
    y = y_full(find(date_quarterly == t1):find(date_quarterly == t2),data_index);
else
    y = dy_full(find(date_quarterly == t1):find(date_quarterly == t2),data_index);
    for j = 1:length(differenced_index_full)
        differenced_index = [differenced_index find(data_index == differenced_index_full(j))];
    end
end

% cut the sample
dateq = date_quarterly(find(date_quarterly == t1):find(date_quarterly == t2));
y_nm = y_nm_full(data_index);

% demean the labor variable
if ismember("labor",y_nm)
    y(:,find(y_nm == "labor")) = y(:,find(y_nm == "labor")) - mean(y(:,find(y_nm == "labor")));
end

% plot the sample
if 0
    for j = 1:length(data_index)
        subplot(3,3,j)
        plot(dateq,y(:,j),'.-','linewidth',1.2)
        recessionplot
        axis tight
        title(y_nm(j))
    end
end

% some steady state values

% government spending
g = GCI(find(date_quarterly_full == t1):find(date_quarterly_full == t2));

% lump-sum transfer
z = GTRSF(find(date_quarterly_full == t1):find(date_quarterly_full == t2));

% spending to GDP ratio
gy = GCI(find(date_quarterly_full == t1):find(date_quarterly_full == t2))./GDP(find(date_quarterly_full == t1):find(date_quarterly_full == t2));

% transfer to GDP ratio
zy = GTRSF(find(date_quarterly_full == t1):find(date_quarterly_full == t2))./GDP(find(date_quarterly_full == t1):find(date_quarterly_full == t2));

% debt to GDP ratio
by = MKVD(find(date_quarterly_full == t1):find(date_quarterly_full == t2))./GDP(find(date_quarterly_full == t1):find(date_quarterly_full == t2));

fprintf('mean of govt spending share to GDP = %.4f, ',mean(gy))
fprintf('mean of market value of debt to GDP = %.4f\n',mean(by))

% compute the autoregressions
lz = log(zy);
lx = lz(1:end-1);
ly = lz(2:end);
a = lx'*ly/(lx'*lx);
u = ly - a*lx;
ux = u(1:end-1);
uy = u(2:end);
ua = ux'*uy/(ux'*ux);
fprintf('Transfer AR(1) coefficient = %.4f, ', a)
fprintf('Autocorrelation of Transfer AR(1) = %.4f\n', ua)
%% compare primary surplus series
if 1    
    mycol = {'#00539a' '#ff832b' '#198038' '#a56eff'};
    linsty = {'-' '-.' ':' '--'};
    mksty = {'none' 'none' 'none' 'none'};
    %mksty = {'o' 's' 'd' 'x'};
    lw = 3.5;
    figure('name','primary surplus data comparison','WindowState','maximized','color','w')
    tiledlayout(1,1, 'Padding', 'none', 'TileSpacing', 'compact');
    nexttile
    hold on
    p1 = plot(date_quarterly,lnps_identity_NSA,'linewidth',lw,'color',mycol{1},'linestyle',linsty{1},'marker',mksty{1});
    p2 = plot(date_quarterly,lnps_account,'linewidth',lw,'color',mycol{2},'linestyle',linsty{2},'marker',mksty{2});    
    p3 = plot(date_quarterly,-UMP,'linewidth',lw,'color',mycol{3},'linestyle',linsty{3},'marker',mksty{3});
    yline(0,'color','#a2191f');
    hold off
    recessionplot
    set(gca, 'YGrid', 'on', 'XGrid', 'off')
    ticklabel = get(gca,'TickLabel');
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',24,'FontWeight','bold');
    legend([p1 p2 p3],'Identity PS', 'Accounting PS', '-Unemployment','location','south','orientation','horizontal','FontSize',24,'FontName','Consolas','FontWeight','bold');
    legend boxoff
    xlabel('year')
    ylabel('%')
end
%% plot the observable
ftsz = 24;
ftsz2 = 20;
T1 = find(date_quarterly == t1);
T2 = find(date_quarterly == t2);
y = y_full(:,data_index);
y_nm = y_nm_full(data_index);
figure('name','Data for Estimation','WindowState','maximized','color','w')
if length(data_index) <= 9
    tiledlayout(3,3, 'Padding', 'none', 'TileSpacing', 'compact');
else
    tiledlayout(2,7, 'Padding', 'none', 'TileSpacing', 'compact');
end
for j = 1:length(data_index)
    nexttile
    hold on
    plot(date_quarterly(T1:T2),y(T1:T2,j),'linewidth',2,'color',mycol{1})
    recessionplot
    hold off
    axis tight
    ticklabel = get(gca,'TickLabel');
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
    title(y_nm(j),'FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
end

dy = dy_full(:,data_index);
figure('name','Data for Estimation: Differenced','WindowState','maximized','color','w')
if length(data_index) <= 9
    tiledlayout(3,3, 'Padding', 'none', 'TileSpacing', 'compact');
else
    tiledlayout(2,7, 'Padding', 'none', 'TileSpacing', 'compact');
end
for j = 1:length(data_index)
    nexttile
    hold on
    plot(date_quarterly(T1:T2),dy(T1:T2,j),'linewidth',2,'color',mycol{1})
    recessionplot
    hold off
    axis tight
    ticklabel = get(gca,'TickLabel');
    set(gca,'TickLabel',ticklabel,'FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
    title(y_nm(j),'FontName','Consolas','fontsize',ftsz,'FontWeight','bold');
end