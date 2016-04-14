function waveletFigures( T, WLT )
%WAVELETFIGURES Makes some basic analyses & checks using wavelets
% See Achard 2006
% Also Bullmore 2003/04 for Hurst
%
%   Inputs: T,      time series vector (time points / parcels)
%           WLT,    wavelet correlation matrices (cell)
%
% Michael Hart, University of Cambridge, February 2016

%% 1. Time series wavelet breakdowns

% Make wavelet decompositions

WJtX = modwt(T); %la8, conservative levels (both default)

nlevels = size(WJtX, 2);
npoints = size(WJtX, 1);

% Single parcel time series

ts = zeros(npoints,nlevels);
for ilevel = 1:nlevels
     ts(:,ilevel) = WJtX(:,ilevel,20);
end

% Plot
figure1 = figure('Name','Wavelet Scales');
subplot1 = subplot(5,1,1,'Parent',figure1);
hold(subplot1,'on');
plot(ts(:,1),'Parent',subplot1);
title({'Scale 1'});
xlim(subplot1,[0 500]);
set(subplot1,'XTickLabel',{'','','','','','','','','','',''});

subplot2 = subplot(5,1,2,'Parent',figure1);
hold(subplot2,'on');
plot(ts(:,2),'Parent',subplot2);
title('Scale 2');
xlim(subplot2,[0 500]);
set(subplot2,'XTickLabel',{'','','','','','','','','','',''});

subplot3 = subplot(5,1,3,'Parent',figure1);
hold(subplot3,'on');
plot(ts(:,3),'Parent',subplot3);
title('Scale 3');
xlim(subplot3,[0 500]);
set(subplot3,'XTickLabel',{'','','','','','','','','','',''});

subplot4 = subplot(5,1,4,'Parent',figure1);
hold(subplot4,'on');
plot(ts(:,4),'Parent',subplot4);
title('Scale 4');
xlim(subplot4,[0 500]);
set(subplot4,'XTick',[50 100 150 200 250 300 350 400 450 500],'XTickLabel',...
    {'','','','','','','','','',''});

subplot5 = subplot(5,1,5,'Parent',figure1);
hold(subplot5,'on');
plot(ts(:,5),'Parent',subplot5);
xlabel('Time (milliseconds)');
title('Scale 5');
xlim(subplot5,[0 500]);

%% 2. Plot matrices over different wavelet scales

figure2 = figure('Name','Wavelet Scale Matrices');
subplot1 = subplot(2,3,1,'Parent',figure2);
hold(subplot1,'on');
image(WLT{1},'Parent',subplot1,'CDataMapping', 'scaled');
title({'Scale 1'});
xlim(subplot1,[0 116]);
ylim(subplot1,[0 116]);

subplot2 = subplot(2,3,2,'Parent',figure2);
hold(subplot2,'on');
image(WLT{2},'Parent',subplot2,'CDataMapping', 'scaled');
title({'Scale 2'});
xlim(subplot2,[0 116]);
ylim(subplot2,[0 116]);

subplot3 = subplot(2,3,3,'Parent',figure2);
hold(subplot3,'on');
image(WLT{3},'Parent',subplot3,'CDataMapping', 'scaled');
title({'Scale 3'});
xlim(subplot3,[0 116]);
ylim(subplot3,[0 116]);

subplot4 = subplot(2,3,4,'Parent',figure2);
hold(subplot4,'on');
imagesc(WLT{4},'Parent',subplot4, 'CDataMapping', 'scaled');
title({'Scale 4'});
xlim(subplot4,[0 116]);
ylim(subplot4,[0 116]);

subplot5 = subplot(2,3,5,'Parent',figure2);
hold(subplot5,'on');
imagesc(WLT{5},'Parent',subplot5, 'CDataMapping', 'scaled');
title({'Scale 5'});
xlim(subplot5,[0 116]);
ylim(subplot5,[0 116]);


%% 3. Power spectra

N = size(T,2); %number of nodes
ts_spectra=[];
for ii=1:N
    blah=abs(fft(nets_demean(T(:,ii)))); 
    ts_spectra(:,ii)=blah(1:round(size(blah,1)/2));
end

%ts_spectra=mean(ts_spectra,3); %number of subjects, if used

figure('Position',[10 10 1600 900]);

grot=ts_spectra ./ repmat(max(ts_spectra),size(ts_spectra,1),1);
grotm=median(grot,2); %mean of all time points

II=3;  I=ceil(N/II);  gap=0.03;  xw=(1-gap*(II+2))/(II+2);  yh=1-(2*gap);  iii=1;
for i=1:II
  subplot('position', [ gap*i+xw*(i-1) gap xw yh*min( (N-iii+1)/I, 1) ]);
  splurghy=repmat(iii:min(iii+I-1,N),size(grot,1),1);
  plot(splurghy,'color',[0.8 0.8 0.8]); hold on;
  plot(splurghy+repmat(grotm,1,size(splurghy,2)),'color',[0.6 0.6 0.6]); hold on;
  clear grotx;
  for ii=1:I
    if iii<=N
      grotx(:,ii)=grot(:,iii)+iii;
    end
    iii=iii+1;
  end
  plot(grotx); % grid on;
end

subplot('position', [ gap*(II+1)+xw*II gap 2*xw yh ]);
plot(grot); hold on;
plot(grotm,'k','LineWidth',3);

%% 4. Hurst exponent

scaleVariance = zeros(nlevels,1);
for ilevel = 1:nlevels
    scaleVariance(ilevel,1) = var(ts(:,ilevel)');
end

scales = [1 2 3 4 5];
figure; semilogy(scales, scaleVariance, 'o'); 

% H = scaleVariance \ scales;

% or, do things manually
% Calculate line slope
% F = [x1.^0 x1];           % make design matrix [1,x]
% c = F\y1;                 % get least-squares fit
% res = y1 - F*c;           % calculate residuals
% r2 = 1 - var(res)/var(y); % calculate R^2


end

