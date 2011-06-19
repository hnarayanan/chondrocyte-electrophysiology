%%
%% Estimate parameters for the A->B->C model
%% jbr,  11/2008

clear all
close all

% ensure sundials/cvode is in the path
if (~ exist ('CVode'))
  if (exist ('startup_STB'))
    startup_STB;
  else
    error ('%s: requires the Sundials Toolbox, please see the README file for installation instructions');
  end
end

global model 

k1 = 2;
k2 = 1;
ca0 = 1;
cb0 = 0;
cc0 = 0;
thetaac = [k1; k2];


tfinal = 6;
nplot = 100;
model.odefcn  = @massbal_ode;
model.tplot   = linspace(0, tfinal, nplot)';
model.param   = thetaac;
model.ic = [ca0; cb0; cc0];

objective.estflag = [1, 2];
objective.paric   = [0.5; 3];
objective.parlb   = [1e-4; 1e-4;];
objective.parub   = [10; 10];

measure.states = [1,2,3];
%% load measurements from a file
table = load ('ABC_data.dat')
measure.time = table(:,1);
measure.data = table(:,1+measure.states);

%% estimate the parameters
estimates = parest(model, measure, objective);


disp('Estimated Parameters and Bounding Box')
[estimates.parest estimates.bbox]


%%plot the model fit to the noisy measurements
figure(1);
plot(model.tplot, estimates.x, measure.time, measure.data, 'o');

%% if 2-d, plot confidence ellipse
if (length(objective.estflag) == 2)
  [xe, ye] = ellipse(estimates.Pinv, estimates.b, 100, estimates.parest);
  figure(2); 
  plot(xe, ye)
endif

table = [model.tplot, estimates.x'];
save ABC.dat table
