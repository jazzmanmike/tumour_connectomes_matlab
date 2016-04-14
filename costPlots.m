%costplots script

figure1 = figure('Name','Measure cost plots');

subplot1 = subplot(5,2,1,'Parent',figure1);
hold(subplot1,'on');
plot(squeeze(mean(nodalMetrics(:,1,:))),'Parent',subplot1);
title({'Degree'});

subplot2 = subplot(5,2,2,'Parent',figure1);
hold(subplot2,'on');
plot(globalMetrics(:,1),'Parent',subplot2);
title({'Strength'});

subplot3 = subplot(5,2,3,'Parent',figure1);
hold(subplot3,'on');
plot(globalMetrics(:,4),'Parent',subplot3);
title({'Clustering'});

subplot4 = subplot(5,2,4,'Parent',figure1);
hold(subplot4,'on');
plot(squeeze(mean(nodalMetrics(:,4,:))),'Parent',subplot4);
title({'Closeness'});

subplot5 = subplot(5,2,5,'Parent',figure1);
hold(subplot5,'on');
plot(globalMetrics(:,5),'Parent',subplot5);
title({'Global efficiency'});

subplot6 = subplot(5,2,6,'Parent',figure1);
hold(subplot6,'on');
plot(globalMetrics(:,6),'Parent',subplot6);
title({'Betweenness'});

subplot7 = subplot(5,2,7,'Parent',figure1);
hold(subplot7,'on');
plot(squeeze(mean(nodalMetrics(:,6,:))),'Parent',subplot7);
title({'Z score'});

subplot8 = subplot(5,2,8,'Parent',figure1);
hold(subplot8,'on');
plot(squeeze(mean(nodalMetrics(:,7,:))),'Parent',subplot8);
title({'Participation'});

subplot9 = subplot(5,2,9,'Parent',figure1);
hold(subplot9,'on');
plot(squeeze(mean(nodalMetrics(:,8,:))),'Parent',subplot9);
title({'Eigenvector'});

subplot10 = subplot(5,2,10,'Parent',figure1);
hold(subplot10,'on');
plot(globalMetrics(:,7),'Parent',subplot10);
title({'Semi-metricity'});
