function metricsSides( nodalMetrics, side )
%METRICSSIDES Plots relationships between sides of measures
%   
%   metricsSides(nodalMetrics, side)
%
%   Inputs: nodalMetrics,   array of metrics e.g. from myHeavyMeasures
%           side,           vector of side indices
%
% Michael Hart, University of Cambridge, February 2017

%% Define 

otherSide = setdiff(1:116, side);

%% Plot 

figure1 = figure('Name','Left versus Right measure comparisons');

subplot1 = subplot(3,4,1,'Parent',figure1);
hold(subplot1,'on');
plot(mean(nodalMetrics(side,:,1),2), mean(nodalMetrics(otherSide,:,1),2),'o');
title({'Degree'});
xlim(subplot1,[0 120]);
ylim(subplot1,[0 120]);
%R = corr(mean(nodalMetrics(left,:,1),2), mean(nodalMetrics(right,:,1),2));

subplot2 = subplot(3,4,2,'Parent',figure1);
hold(subplot2,'on');
plot(mean(nodalMetrics(side,:,2),2), mean(nodalMetrics(otherSide,:,2),2),'o');
title({'Strength'});
ylim(subplot2,[0 80]);

subplot3 = subplot(3,4,3,'Parent',figure1);
hold(subplot3,'on');
plot(mean(nodalMetrics(side,:,3),2), mean(nodalMetrics(otherSide,:,3),2),'o');
title({'Clustering'});
ylim(subplot3,[0.2 0.6])

subplot4 = subplot(3,4,4,'Parent',figure1);
hold(subplot4,'on');
plot(mean(nodalMetrics(side,:,4),2), mean(nodalMetrics(otherSide,:,4),2),'o');
title({'Normalised clustering'});
ylim(subplot4,[0.2 0.6])

subplot5 = subplot(3,4,5,'Parent',figure1);
hold(subplot5,'on');
plot(mean(nodalMetrics(side,:,5),2), mean(nodalMetrics(otherSide,:,5),2),'o');
title({'Local efficiency'});
ylim(subplot5,[0.2 0.6])

subplot6 = subplot(3,4,6,'Parent',figure1);
hold(subplot6,'on');
plot(mean(nodalMetrics(side,:,6),2), mean(nodalMetrics(otherSide,:,6),2),'o');
title({'Closeness'});
ylim(subplot6,[0.2 0.6])

subplot7 = subplot(3,4,7,'Parent',figure1);
hold(subplot7,'on');
plot(mean(nodalMetrics(side,:,7),2), mean(nodalMetrics(otherSide,:,7),2),'o');
title({'Betweenness'});
xlim(subplot7,[0 160]);
ylim(subplot7,[0 160]);

subplot8 = subplot(3,4,8,'Parent',figure1);
hold(subplot8,'on');
plot(mean(nodalMetrics(side,:,8),2), mean(nodalMetrics(otherSide,:,8),2),'o');
title({'Z-score'});

subplot9 = subplot(3,4,9,'Parent',figure1);
hold(subplot9,'on');
plot(mean(nodalMetrics(side,:,9),2), mean(nodalMetrics(otherSide,:,9),2),'o');
title({'Participation'});
ylim(subplot9,[0.3 0.7]);

subplot10 = subplot(3,4,10,'Parent',figure1);
hold(subplot10,'on');
plot(mean(nodalMetrics(side,:,10),2), mean(nodalMetrics(otherSide,:,10),2),'o');
title({'Eigenvector'});
ylim(subplot10,[0 0.15]);

subplot11 = subplot(3,4,11,'Parent',figure1);
hold(subplot11,'on');
plot(mean(nodalMetrics(side,:,11),2), mean(nodalMetrics(otherSide,:,11),2),'o');
title({'Pagerank'});

subplot12 = subplot(3,4,12,'Parent',figure1);
hold(subplot12,'on');
plot(mean(nodalMetrics(side,:,12),2), mean(nodalMetrics(otherSide,:,12),2),'o');
title({'Semi-metricity'});
ylim(subplot12,[0.6 1]);

end

