% SurrModelEDFNonnormal
% This example creates some slightly troubled data and shows that the
% problems propagate beyond the troubled region.  This demonstrates that
% fitting an EDF to converge to a CDF is difficult because the errors do
% not have the desired distribution.
% For this example we reuse the Generalized Pareto from
% SurrModelEDF1D

% Standardize the random results
global GAUSSQR_PARAMETERS
GAUSSQR_PARAMETERS.RANDOM_SEED(0);

% Create some slightly off data, with an extra 10% of the points less
% than the median
x = sort([rand(20,50)*(2-sqrt(2));...
          icdf('gp',rand(200,50),-.5,1,0)]);

% Form the EDF values at the data points
EDFs = repmat((1:220)'/220,1,50);

% Compute the difference between the EDF and CDF
% Only store the results for x>1 (outside of bad region)
ECdiff = EDFs(x>1)-cdf('gp',x(x>1),-.5,1,0);
n = length(ECdiff);

% Plot the results in a qqplot
h_fig = figure;
h = normplot(ECdiff);
legend(h([1,3]),'Data','Normal Expectation','location','northwest')

% Compute the skewness and add it to the plot
skewness = sqrt(n^2-n)/(n-2)*...
           sum((ECdiff-mean(ECdiff)).^3)/n/std(ECdiff)^3;
annotation('textbox',[.6,.2,.2,.07],'string',...
           sprintf('skewness = %3.2f',skewness'))
       
% Run the KS test to see if it fits a normal distribution
[h,p] = kstest((ECdiff-mean(ECdiff))/std(ECdiff))