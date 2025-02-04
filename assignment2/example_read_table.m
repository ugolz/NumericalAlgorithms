
TT=readtable('covid19_regional_data.xlsx','Sheet','Dead'); % TT is a table structure
data_dead=TT.Marche;     % y is a vector
t=datenum(TT.date); % transform dates in a numerical format
clear TT

TT=readtable('covid19_regional_data.xlsx','Sheet','NewPositive'); % TT is a table structure
data_positive=TT.Marche;     % y is a vector
clear TT

TT=readtable('covid19_regional_data.xlsx','Sheet','Recovered'); % TT is a table structure
data_recovered=TT.Marche;     % y is a vector
clear TT

%% plot
figure()
set(gcf,'color','w')
plot(t,log(data_dead),'o'),
legend('data','Location','northwest');

datetick('x','mmm-dd')
title('Marche')
ylabel('new positive cases')           