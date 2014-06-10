clc
close all


plot(1,1,'r.', 2,2', 'b.', 3,3, 'k.');
legendHandle = legend('plot 1', 'plot 2 with longer title', 'plot 3');

tableHead = {'\theta_0' '\phi' 'df/dx'};
tableContent = rand(3);

legendTable(legendHandle, tableHead, tableContent);


%%

clc

legendlocation = 'NorthWestOutside'

regexprep(legendlocation, 'inside|outside$', '', 'ignorecase')