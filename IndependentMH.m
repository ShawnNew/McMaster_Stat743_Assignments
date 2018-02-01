%------------------------------------------------------
% EXAMPLE 1: METROPOLIS-HASTINGS
% BLOCK-WISE SAMPLER (BIVARIATE NORMAL)
rand('seed' ,12345);
 
D = 2; % # VARIABLES
nBurnIn = 100;
 
% TARGET DISTRIBUTION IS A 2D NORMAL WITH STRONG COVARIANCE
p = inline('mvnpdf(x,[0 0],[1 0.8;0.8 1])','x');
 
% PROPOSAL DISTRIBUTION STANDARD 2D GUASSIAN
q = inline('mvnpdf(x,mu)','x','mu')
 
nSamples = 5000;
minn = [-3 -3];
maxx = [3 3];
 
% INITIALIZE BLOCK-WISE SAMPLER
t = 1;
x = zeros(nSamples,2);
x(1,:) = randn(1,D);
 
% RUN SAMPLER
while t < nSamples
    t = t + 1;
 
    % SAMPLE FROM PROPOSAL
    xStar = mvnrnd(x(t-1,:),eye(D));
 
    % CORRECTION FACTOR (SHOULD EQUAL 1)
    c = q(x(t-1,:),xStar)/q(xStar,x(t-1,:));
 
    % CALCULATE THE M-H ACCEPTANCE PROBABILITY
    alpha = min([1, p(xStar)/p(x(t-1,:))]);
 
    % ACCEPT OR REJECT?
    u = rand;
    if u < alpha
        x(t,:) = xStar;
    else
        x(t,:) = x(t-1,:);
    end
end
 
% DISPLAY
nBins = 20;
bins1 = linspace(minn(1), maxx(1), nBins);
bins2 = linspace(minn(2), maxx(2), nBins);
 
% DISPLAY SAMPLED DISTRIBUTION
ax = subplot(121);
bins1 = linspace(minn(1), maxx(1), nBins);
bins2 = linspace(minn(2), maxx(2), nBins);
sampX = hist3(x, 'Edges', {bins1, bins2});
hist3(x, 'Edges', {bins1, bins2});
view(-15,40)
 
% COLOR HISTOGRAM BARS ACCORDING TO HEIGHT
colormap hot
set(gcf,'renderer','opengl');
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
xlabel('x_1'); ylabel('x_2'); zlabel('Frequency');
axis square
set(ax,'xTick',[minn(1),0,maxx(1)]);
set(ax,'yTick',[minn(2),0,maxx(2)]);
title('Sampled Distribution');
 
% DISPLAY ANALYTIC DENSITY
ax = subplot(122);
[x1 ,x2] = meshgrid(bins1,bins2);
probX = p([x1(:), x2(:)]);
probX = reshape(probX ,nBins, nBins);
surf(probX); axis xy
view(-15,40)
xlabel('x_1'); ylabel('x_2'); zlabel('p({\bfx})');
colormap hot
axis square
set(ax,'xTick',[1,round(nBins/2),nBins]);
set(ax,'xTickLabel',[minn(1),0,maxx(1)]);
set(ax,'yTick',[1,round(nBins/2),nBins]);
set(ax,'yTickLabel',[minn(2),0,maxx(2)]);
title('Analytic Distribution')