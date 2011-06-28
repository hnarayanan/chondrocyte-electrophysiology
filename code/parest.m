## Copyright (C) 2007 John W. Eaton
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, write to the Free
## Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.

function retval = parest(model, measure, objective)
%  Function File:	estimates = @parest(model, measure, objective)
%  
%  Estimates parameters p for differential equation models
%		dx/dt = f(x,t,p)
%		x(0) = x0(p)
%  for given masurements y_i = h(x(t_i)) at time points t_i.
% 
%  The estimates are obtained by solving a weighted least squares problem
% 
%		min Sum ( h(x(t_i)) - y_i )' W ( h(x(t_i)) - y_i )
%		 p   i
%
%  with weight matrix W.
%
%  Input Arguments
%  ===============
%  The three input arguments model, measure, objective 
%  are fields and contain:
%
%  1. model
%  --------
%  model.odefcn (required)
%		Provides the right hand side of the 
%		differential equation in the form 
%		xdot = f(x, t, parameters).
%
%  model.cvodesfcn (optional)
%		Provides the right hand side of the 
%		differential equation in the form 
%		expected by CVODES,
%		[xdot, flag,new_data] = cvodesrhs(t, x, data).
%		If you provide this function and model.odefcn,
%		then the CVODES RHS function is used and 
%		model.odefcn is not called.
%
%  model.dodedp, model.dodedx (recommended)
%  		Partial derivatives of f with respect to the 
%		model parameters and the state respectively.
%  		Necessary for exact calculation of gradients 
%		of the objective function.
%  
%  model.ic (required)
%  		Column vector with initial condition of the 
%		system states.
%		If several data sets for different experiments 
%		corresponding to different initial conditions 
%		are given, this is a matrix in which each the
%		i-th column is the state initial condition 
%		corresponding to the i-th experiment.
%
%  model.param (required)
%  		Values of model parameters or initial guesses for 
%		the parameters which are estimated.
%  
%  model.stateic_par (optional)
%  		Vector or matrix of the same size as model.ic
%		which describes which initial conditions of
%		the state are parameters (and possibly have 
%		to be estimated).
%		E.g. stateic_par = [0; 1; 0; 2]
%		states that the second initial condition is
%		parameter 1 and the fourth initial condition
%		is parameter 2.
%  
%  model.tplot (optional)
%  		Time points for which estimated states are 
%		returned, default: measure.time (cf. below).
%
%  model.atol
%  model.rtol (optional)
%		Absolute and relative tolerances for the
%		integration of the differential equations
%		with CVodes. Default: sqrt(eps)
%  model.odesteps (optional)
%		Maximum number of steps for solving the 
%		differential equations with CVodes.
%		Default: 100,000
%
%  2. measure
%  ----------
%  measure.states
%  or
%  measure.statefcn, measure.dstatedx
%  or
%  measure.measurematrix (recommended)
%		The user can either provide the measured states
%		as a vector or a function y = h(x) along with
%		its derivative dh/dx or the measurematrix C
%		for y = C*x.
%		If neither is provided, parest assumes all 
%		states to be measured.
%		If two or more are provided, an error is given.
%
%  measure.time (required)
%		Column vector of the measurement times.
%		If several data sets for different initial 
%		conditions are given, this is a matrix with
%		the second dimension corresponding to the 
%		index of the experiment.
%  
%  measure.data (required)
%  		Matrix of data, the rows correspond to the 
%		measured states or output functions h 
%		respectively and the columns correspond to 
%		the times provided in measure.time.
%		If several data sets for different initial 
%		conditions are given, this is a 3-dimensional 
%		structure with the third dimension corresponding 
%		to the index of the experiment.
%
%  measure.weight (optional)
%		Matrix W for weighted least squares or covariance
%		matrix of measurement noise, default: eye(n).
%  
%  measure.alpha (optional)
% 		Alpha for confidence intervals, default: 0.95
%  
%  3. objective
%  ------------
%  objective.estflag (recommended)
%  		Vector of indices of parameters in model.param 
%		that have to be estimated.
%		If estflag is not given, all parameters are 
%		estimated. If estflag is an empty vector [], 
%		the output values are calculated for the given 
%		initial guess of parameters (cf. objective.paric).
%  
%  objective.paric, objective.parlb, objective.parub (recommended)
%  		Vectors of length of estflag with initial guess, 
%		lower and upper bounds for the parameters which are 
%		estimated. When no estflag is provided, these vectors 
%		have to be of the same length as model.param.
%		For empty / missing vectors the following defaults 
%		are used:
%		Default for objective.paric: values in model.param
%		Default for objective.parlb: realmin
%		Default for objective.parub: realmax
%		If only scalars are given for objective.parlb and
%		objective.parub, those bounds are used for all 
%		parameters.
%
%  objective.logtr (optional)
%  		Vector of logicals of length of estflag which 
%		indicates for which parameters a logarithmic tranformation 
%		should be performed.
%  
%  objective.sensitivityAdjointMethod (optional)
%  		Logical, if true use adjoint method for gradient
%		calculations, if false forward sensitivity 
%		analysis is used.
%
%  objective.printlevel (optional)
%		Decides which informations about the current
%		optimization process are displayed on the screen.
%		The options are currently available:
%		0: (default) no output on screen is given
%		1: parameter values, value of objective function and
%			gradient are displayed on screen during
%			optimization
%
%  Return Values
%  =============
%  The return value estimates is a field and contains:
%
%  estimates
%  ---------
%  estimates.parest, estimates.parestlogtr
%  		Estimated parameter values or its logarithm for
%  		those specified in objective.logtr.
%  
%  estimates.obj
%  		Value of the objective function.
%  
%  estimates.grad, estimates.gradlogtr
%  estimates.Hgn, estimates.Hgnlogtr
%  estimates.HFD
%		Gradient and Hessian of objective function with 
%		respect to the parameters or their logarithms for
%		those specified in objective.logtr.
%		The Hessian in Hgn and Hgnlogtr is calculated by 
%		using a Gauss-Newton Approximation, whereas HFD 
%		is calculated with a two sided finite difference 
%		approach applied on the gradients.
%
%  estimates.info  
%		Information flag returned by the optimization program.
%  
%  estimates.mu
%		Matrix of adjoint variable in form mu(:,measure.time)
%  		if adjoint method is used for gradient calculation
%		(cf. objective.sensitivityAdjointMethod).
%		If several data sets for different initial conditions 
%		are given, this is a 3-dimensional structure with the
%		third dimension corresponding to the index of the 
%		experiment.
%  
%  estimates.x
%  estimates.sx
%  		State and sensitivities for the optimal parameter
%		values at time tplot or by default the measurement
%		times.
%		x is a matrix with the ith column being the state at
%		the ith time point, i.e.
%		x(:,i) = x(t_i).
%		sx is a three-dimensional structure, with the last
%		argument corresponding to time, i.e.
%		sx(:,:,i) = dx/dp(t_i).
%		If several data sets for different initial conditions
%		are given, x (sx) is a 3 (4)-dimensional structure
%		with the third (fourth) dimension corresponding to
%		the index of the experiment.
%
%  estimates.y
%		Output of system for the optimal parameter
%		values at time tplot or by default the measurement
%		times (cf. estimates.x).
%
%  estimates.measpred
%  		Matrix of predicted measurements with columns
%		corresponding to the measurement times.
%  estimates.resid
%  		Matrix of residuals in the same format as 
%		estimates.measpred.
%		If several data sets for different initial 
%		conditions are given, resid and measpred are 
%		3-dimensional structures with the third dimension 
%		corresponding to the index of the experiment.
%
%  estimates.measvar
%  		Estimate of measurement error covariance matrix.
%
%
%  estimates.bbox, estimates.bboxFD
%		Bounding box for confidence intervals for given 
%		measure.alpha confidence level using the two methods of 
%		finding the Hessian matrix.
%
%  estimates.Pinv, estimates.b
%               Confidence ellipse:
%  (theta - estimates.parest)' * (estimates.Pinv) * (theta - estimates.parest)  \leq b
%
%  See also: -

  global mod_ meas_ obj_
  global solveroption

  %% create local copies of these objects to pass to likelihood
  meas_  = measure; obj_ = objective;  mod_ = model;

  %% 1. Check model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if(~isfield(mod_,'odefcn') && ~isfield(mod_,'cvodesfcn'))
	error('parest: model.odefcn and model.cvodesfcn missing: no rhs of ode provided');
  end

  if(~isfield(mod_,'param'))
	error('parest: model.param missing');
  end

  if(~isfield(mod_,'ic'))
	error('parest: model.ic missing');
  end

  %% if user provides no stateic_par or par_stateic
  if(~isfield(mod_,'stateic_par') || isempty(mod_.stateic_par))
	mod_.stateic_par = zeros(size(mod_.ic));
  end

  %% if entry in model.stateic_par is too large
  if(max(mod_.stateic_par) > numel(mod_.param))
	error('parest: entry in model.stateic_par exceeds number of parameters')
  end

  %% if sizes don't match
  %% e.g. if only model.stateic_par(3) = 2 was used
  %% for a system with 4 states
  %% -> fill up with zeros
  if (~isequal( size(mod_.ic), size(mod_.stateic_par) ) )
	warning('parest: sizes of model.ic and model.stateic_par do not match')
	stateic_par = zeros(size(mod_.ic));
	for i = 1:size(stateic_par,1)
		for j = 1:size(stateic_par,2)
			stateic_par(i,j) = mod_.stateic_par(i,j);
		end
	end
	mod_.stateic_par = stateic_par;
  end


  if(isfield(mod_,'tplot'))
  if(size(mod_.tplot,1) == 1 && size(mod_.tplot,2) ~= 1 )
	% if tplot is provided as row vector -> transpose it
	mod_.tplot = mod_.tplot';
	warning('parest: model.tplot is row vector')
  end
  end

  if(~isfield(mod_,'atol'))
	% absolute tolerance for CVodes
	mod_.atol = sqrt(eps);
  elseif(~isscalar(mod_.atol))
	% absolute tolerance for CVodes
	error('parest: model.atol has to be scalar')
  end

  if(~isfield(mod_,'rtol'))
	% relative tolerance for CVodes
	mod_.rtol = sqrt(eps);
  elseif(~isscalar(mod_.rtol))
	% relative tolerance for CVodes
	error('parest: model.rtol has to be scalar')
  end

  if(~isfield(mod_,'odesteps'))
	% maximum number of steps for CVodes
	mod_.odesteps = 100000;
  elseif(~isscalar(mod_.odesteps))
	% max number of steps tolerance for CVodes
	error('parest: model.atol has to be scalar')
  end

  %% 2. Check objective %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% if user provides no estflag, estimate all parameters
  if(~isfield(obj_,'estflag'))
	obj_.estflag = 1:numel(mod_.param);
	warning('parest: objective.estflag is missing. All parameters are estimated')
  end
  
  %% if estflag is empty, don't estimate any parameters 
  %% but calculate objective function etc
  if(isempty(obj_.estflag))
	obj_.estflag = 1:numel(mod_.param);
	estflagempty = true;
	warning('parest: objective.estflag is empty. No optimization is performed')
  else
	estflagempty = false;	
  end

  %% make column vector
	if ( size(obj_.estflag,2) > 1 )
		obj_.estflag = obj_.estflag';
	end

  %% number of parameters to estimate
  nparest = numel(obj_.estflag);

  %% if user provides no special initial guess take values out of model.param
  if(~isfield(obj_,'paric') || isempty(obj_.paric) )
	obj_.paric = mod_.param(obj_.estflag);
	warning('parest: objective.paric is missing. Model.param used as initial guess')
  end

  %% if user doesn't provide logtr
  if(~isfield(obj_,'logtr') || isempty(obj_.logtr))
	%as default don't use log transform
	obj_.logtr = logical(zeros(nparest,1));
  end
		
  %% check if sizes in objective are correct
  if (numel(obj_.logtr) ~= nparest)
	error('parest: Number of elements in objective.logtr must be equal to length of objective.estflag') 
  end
  if (~isfield(obj_,'parlb') || isempty(obj_.parlb))
	% if no lower bound is given
	obj_.parlb = realmin*ones(nparest,1);%1e-5*obj_.paric;
	warning('parest: objective.parlb is missing')
	%% error('parest: objective.parlb missing')
  elseif (numel(obj_.parlb) ~= nparest)
	if (isscalar(obj_.parlb))
		obj_.parlb = obj_.parlb*ones(1,nparest);
	else
	error('parest: Number of elements in objective.parlb must be equal to length of objective.estflag')
	end
  end
  if (~isfield(obj_,'parub') || isempty(obj_.parub))
	% if no upper bound is given
	obj_.parub = realmax*ones(nparest,1);%1e5*obj_.paric;
	warning('parest: objective.parub is missing')
	%% error('parest: objective.parub missing')
  elseif (numel(obj_.parub) ~= nparest)
	if (isscalar(obj_.parub))
		obj_.parub = obj_.parub*ones(1,nparest);
	else
	error('parest: Number of elements in objective.parub must be equal to length of objective.estflag')
	end
  end
  % ensure that parub, parlb, and paric are column vectors
  obj_.parub = obj_.parub(:);
  obj_.parlb = obj_.parlb(:);
  obj_.paric = obj_.paric(:);
  % ensure ub > lb
  if (any( obj_.parub - obj_.parlb <= 0))
	error('parest: Upper bound has to be greater than lower bound')
  end
  if (numel(obj_.paric) ~= nparest)
	error('parest: Number of elements in objective.paric must be equal to length of objective.estflag')
  end

  if(~isfield(obj_,'printlevel'))
	obj_.printlevel = 0;
  end

  %% 3. Check measure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if(isfield(meas_,'states') + isfield(meas_,'statefcn') + isfield(meas_,'measurematrix') > 1)
	error('parest: more than one of measure.statefcn, states and measurematrix provided')
  end

  if(~isfield(meas_,'states') && ~isfield(meas_,'statefcn') && ~isfield(meas_,'measurematrix'))
	meas_.states = 1:size(mod_.ic,1);
	warning('parest: none of measure.statefcn, states and measurematrix provided')
  end

  if(isfield(meas_,'statefcn') && ~isfield(meas_,'dstatedx'))
	error('parest: measure.dstatedx missing');
  end

  if(~isfield(meas_,'data'))
	error('parest: measure.data missing');
  end

  if(~isfield(meas_,'time'))
	error('parest: measure.time missing');
  end

  if(size(meas_.time,1) == size(mod_.ic,2) && size(meas_.time,2) ~= size(mod_.ic,2))
	% if measurement time is provided as column vector -> transpose it
	meas_.time = meas_.time(:);
  end

  %% if transpose of measurement.data is given
  %% data is supposed to be in format (states,time)
  if (size(meas_.data,2) ~= size(meas_.time,1) && size(meas_.data,1) == size(meas_.time,1)) 
	for k = 1:size(meas_.data,3)
	measdata(:,:,k) = (meas_.data(:,:,k))';
	end
	meas_.data = measdata;
	warning('parest: measure.data seems to be transposed')
  end

  if (size(meas_.time,2) ~= size(meas_.data,3) || size(meas_.time,2) ~= size(mod_.ic,2) || size(meas_.time,1) ~= size(meas_.data,2))
	error('parest: sizes of measure.data, measure.time and model.ic are not compatible')
  end

  %% number of measurements at each time point, e.g. =2 if x1 and x2 are measured
  nmeas = size(meas_.data, 1);

  if(isfield(meas_,'states'))
	if( numel(meas_.states) ~= nmeas)
		error('parest: sizes of data and measurement.states do not match')
	end
  elseif(isfield(meas_,'statefcn'))
	y0 = meas_.statefcn(mod_.ic);
	dy0dx = meas_.dstatedx(mod_.ic);
	if( size(y0,1) ~= nmeas)
		error('parest: sizes of data and measurement.statefcn do not match')
	end
	if( size(dy0dx,1) ~= nmeas || size(dy0dx,2) ~= numel(mod_.ic) )
		error('parest: size of measurement.dstatedx seems wrong')
	end
  elseif(isfield(meas_,'measurematrix'))
	if( size(meas_.measurematrix,1) ~= nmeas || size(meas_.measurematrix,2) ~= numel(mod_.ic) )
		error('parest: size of measurematrix seems wrong')
	end
  end

  %% ensure that measure weight matrix is symmetric
  %% and has appropriate size
  if(isfield(meas_,'weight'))
	if (size(meas_.weight) == nmeas*ones(1,2))
  		meas_.weight = 0.5*(meas_.weight+meas_.weight');
	else
		error('parest: size of measurement.weight does not fit to measurement.data')
	end
  else
	meas_.weight = eye(nmeas);
	%warning('parest: weight matrix missing. Identity matrix used')
  end
	
  if(size(mod_.ic,2) ~= size(meas_.data,3))
	error('parest: number of data sets with different ICs does not fit with meaurement.data')
  end

  if(size(mod_.ic,2) ~= size(meas_.time,2))
	error('parest: number of initial condition sets does not fit with meaurement.time')
  end

  %% 4. prepare optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% log transform the specified parameters
	% ensure objective.logtr to be logical
	obj_.logtr = logical(obj_.logtr);
	trind =  obj_.logtr;
	obj_.paric(trind) = log(obj_.paric(trind));
	obj_.parlb(trind) = log(obj_.parlb(trind));
	obj_.parub(trind) = log(obj_.parub(trind));

  %% initialize the par sensitivities; note the pars that are ICs
  %% sx0 = 0 for parameters that are not state IC
  %% sx0 = 1 for parameters that are state IC
  %% third dimension corresponds to different experiments
	mod_.sx0 = zeros(size(mod_.ic,1), numel(obj_.estflag), size(mod_.ic,2));
	
	for i = 1:size(mod_.ic,1) %for all states
		for k = 1:size(mod_.ic,2) % for all experiments
			if (mod_.stateic_par(i,k) > 0)
				mod_.sx0(i,mod_.stateic_par(i,k) ,k) = 1;
			end
		end %for k = ...
	end % for i = ...

  %% 5. call optimizer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(estflagempty)
	%% if estflag is empty -> no parameters to estimate, no optimization
	theta = obj_.paric;
	info = 'estflag empty, no optimization performed';
  else
	%% optimize in order to
	%% find best fit of parameters
	%% solveroption determines if/how gradient is provided for optimization
	solveroption = 2;
	if (isfield(obj_,'sensitivityAdjointMethod') && obj_.sensitivityAdjointMethod)
	        if(isfield(mod_,'dodedx') && isfield(mod_,'dodedp'))
			solveroption = 3;
        	else
        		warning('parest: adjoint method cannot be used because model.dodedp or dodedx missing')
        	end
	end

	%% minimize negative logarithm of likelihood function
	switch solveroption
	case 1 % solve without provided gradients - currently not used
	  [theta, phi, info] = fmincon (@likelihood, obj_.paric, ...
        		[], [], [], [], obj_.parlb, obj_.parub);
	case 2 % solve using forward sensitivities
	  options = optimset('GradObj','on');
	  [theta, phi, info] = fmincon ({@likelihood,@gradi}, obj_.paric, ...
        		[], [], [], [], obj_.parlb, obj_.parub,[],options);
	case 3	% solve using adjoint method for gradients
	  options = optimset('GradObj','on');
	  [theta, phi, info] = fmincon ({@likelihoodA,@gradi}, obj_.paric, ...
        		[], [], [], [], obj_.parlb, obj_.parub,[],options);
	end %end switch solveroption
  end %if(estflagempty)

  %% 6. work on return values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %% (Insurance) call to likelihood with optimal parameters to compute solution,
  %% sensitivities, gradient and Hessian. Not necessary if optimizer's
  %% last call was with optimal values.
	phi = likelihood(theta);

  % if tlot is given, solve ODE for all time points in tplot
  if(isfield(mod_,'tplot'))
	data.theta = mod_.param;
	%% set the model parameters; note the log transformation
	data.theta(obj_.estflag) = theta;
	trind =  obj_.estflag(obj_.logtr);
	data.theta(trind) = exp(theta(obj_.logtr));
	%% solve the model and sensitivites at tplot
		if (~isfield(mod_,'cvodesfcn'))
			if (~isfield(mod_,'dodedx') | ~isfield(mod_,'dodedp') )
				[x, sx, status] = cvode_sens(@oderhs, [], mod_.ic, ...
					obj_.estflag, data, mod_.sx0, mod_.tplot);
			else
				[x, sx, status] = cvode_sens(@oderhs, @sensrhs, mod_.ic, ...
					obj_.estflag, data, mod_.sx0, mod_.tplot);  
			end
		else
			if (~isfield(mod_,'dodedx') | ~isfield(mod_,'dodedp') )
				[x, sx, status] = cvode_sens(mod_.cvodesfcn, [], mod_.ic, ...
				obj_.estflag, data, mod_.sx0, mod_.tplot);
			else
				[x, sx, status] = cvode_sens(mod_.cvodesfcn, @sensrhs, mod_.ic, ...
					obj_.estflag, data, mod_.sx0, mod_.tplot);  
			end
		end		
	% if cvode_sens failed
	if (status ~= 0)
		error('parest: CVode failed with status = %d', status);
	end
	retval.x = x;
	retval.sx = sx;
  else
	%% if no tplot is provided by user return estimates and 
	%% sensitivities at available measurement time point
	retval.x = obj_.x;
	retval.sx = obj_.sx;
  end %if(isfield(mod_,'tplot'))

  nics = size(retval.x, 3);
  ntimes = size(retval.x,2);
  for k = 1:nics %for each IC = dataset
	for i = 1: ntimes
	if (isfield(meas_,'states'))
		% if user provides vector of measured states
		output(:,i,k) = retval.x(meas_.states,i,k);
	elseif (isfield(meas_,'statefcn'))
		% if user provides a measurement state function
		output(:,i,k) = meas_.statefcn(retval.x(:,i,k));
	elseif (isfield(meas_,'measurematrix'))
		% if output matrix y = measurematrix*x is provided
		output(:,i,k) = meas_.measurematrix*(retval.x(:,i,k));
	end
	end % for i = ...
  end % for k = ...

  %% store output for return
	retval.y = output;
	retval.obj  = phi;
	retval.info = info;
  %% log transformed parameters
	retval.parestlogtr = theta;
  %% gradient and Hessian; pull these out of common; computed in likelihood
	retval.resid = obj_.resid;
	retval.measpred = obj_.measpred;
	retval.grad  = obj_.grad;
	retval.gradlogtr = obj_.gradlogtr;
  %% undo log transform the specified parameters
	thetawolog = theta;
	thetawolog(obj_.logtr) = exp(theta(obj_.logtr));
	retval.parest = thetawolog;

  %% if user provides no alpha for confidence intervals
  if(~isfield(meas_,'alpha'))
	meas_.alpha = 0.95;
  end

  if(isfield(obj_,'Hgn'))
	retval.Hgn  = obj_.Hgn;
	retval.Hgnlogtr = obj_.Hgnlogtr;
	%%compute confidence intervals
	p     = length(theta);
	ndata = length(meas_.time);
	Fstat = finv(meas_.alpha, p, ndata-p);
	Hgn   = retval.Hgn;
	if (isscalar(Hgn))
	  invHgn = 1./Hgn;
	else
	  invHgn = inv(Hgn);
	end
        retval.bbox = sqrt(2*p/(ndata-p)*Fstat*phi*diag(invHgn));
        samvar = retval.obj/(ndata-p);
        retval.Pinv = retval.Hgn/(2*samvar);
        retval.b = p*Fstat;
  end

  %% adjoint variable mu if computed in likelihoodA
  %% matrix of form mu(:,time)
  if(isfield(obj_,'mu'))
	retval.mu = obj_.mu;
  end

  %% Estimate of Measurement Variance Matrix
	residv = retval.resid;
	%% first find and delete bad data
	residv(find(~isfinite(residv))) = 0;
	%% formula for Variance
	retval.measvar = 0;
	for k = 1:size(mod_.ic,2)
	retval.measvar = retval.measvar + residv(:,:,k)*residv(:,:,k)'/size(residv,2);
	end

   %% approximate Hessian via Finite Difference
	HFD = zeros(nparest,nparest);
	for j =1:nparest
		FDperturbation = 1*sqrt(1e-3)*thetawolog(j);
		thetaFD = thetawolog;
		thetaFD(j) = thetawolog(j)+FDperturbation; %perturb one parameter
		thetaFD2 = thetawolog;
		thetaFD2(j) = thetawolog(j)-FDperturbation; %perturb one parameter
		thetaFD2(obj_.logtr) = log(thetaFD2(obj_.logtr)); %log trans
		likelihood(thetaFD2);
		grad2 = obj_.grad;
		thetaFD(obj_.logtr) = log(thetaFD(obj_.logtr)); %log trans
		likelihood(thetaFD);
		%HFD(:,j) = (obj_.grad - retval.grad)/FDperturbation;
		HFD(:,j) = (obj_.grad - grad2)/(2*FDperturbation);
    	end
	HFD = 0.5*(HFD+HFD'); % make HFD symmetric
	retval.HFD = HFD;
	  %%compute confidence intervals with FD
	p     = length(theta);
	ndata = length(meas_.time);
	Fstat = finv(meas_.alpha, p, ndata-p);
	if (isscalar(HFD))
	  invHFD = 1./HFD;
	else
	  invHFD = inv(HFD);
    	end
        retval.bboxFD = sqrt(2*p/(ndata-p)*Fstat*retval.obj*diag(invHFD));

end %parest.m
%% END OF FUNCTION PAREST.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phi = likelihood(theta)
%% Calculates negative logarithm of likelihood function,
%% this function is minimized by parest.m.
%% Gradients and sensitivities are calculated by forward calculations
%% using cvode_sens.
  global mod_ meas_ obj_ likelihood_last_theta

  % print on screen if printlevel is set accordingly
  if(obj_.printlevel == 1)
	disp('----------------------');
	disp('examined parameter values:');
	disp(theta);
  end

  %% data.theta contains all parameters, not only those which are estimated!
  data.theta = mod_.param;
  %% set the model parameters; note the log transformation
  data.theta(obj_.estflag) = theta;
  trind =  obj_.estflag(obj_.logtr);
  data.theta(trind) = exp(theta(obj_.logtr));

  %% set the state ICs that are parameters
  for i = 1:numel(mod_.stateic_par)
	if (mod_.stateic_par(i) > 0)
		mod_.ic(i) =  data.theta(mod_.stateic_par(i));
	end
  end

  %% solve the model and sensitivites at the measurement times;
  %% predict measurements
	if (~isfield(mod_,'cvodesfcn'))
		if (~isfield(mod_,'dodedx') || ~isfield(mod_,'dodedp') )
			[x, sx, status] = cvode_sens(@oderhs, [], mod_.ic, ...
			obj_.estflag, data, mod_.sx0,  meas_.time);
		else
			[x, sx, status] = cvode_sens(@oderhs, @sensrhs, mod_.ic, ...
			obj_.estflag, data, mod_.sx0,  meas_.time);  
		end
	else
		if (~isfield(mod_,'dodedx') || ~isfield(mod_,'dodedp') )
			[x, sx, status] = cvode_sens(mod_.cvodesfcn, [], mod_.ic, ...
			obj_.estflag, data, mod_.sx0,  meas_.time);
		else
			[x, sx, status] = cvode_sens(mod_.cvodesfcn, @sensrhs, mod_.ic, ...
			obj_.estflag, data, mod_.sx0,  meas_.time);  
		end
	end
  % if cvode_sens failed
	if (status ~= 0)
		phi = realmax;
		return;
	end

  nmeas = size(meas_.data, 1);
  ntimes = size(meas_.data, 2);
  nics = size(mod_.ic, 2);
  npar = numel(theta);
  resid = zeros(nmeas,ntimes,nics);
  measpred = zeros(nmeas,ntimes,nics);
  nx = size (mod_.ic,1);

  %% is a loop best here?
  phi = 0;
  grad = zeros(npar, 1);
  Hgn  = zeros(npar, npar);
  for k = 1:nics %for each IC = dataset
	for i = 1: ntimes
		% for current time point i and data set k:
	
		if (isfield(meas_,'states'))
			% if user provides vector of measured states
			measpred(:,i,k) = x(meas_.states,i,k);
			resid(:,i,k) = meas_.data(:,i,k) - x(meas_.states,i,k);
			IdentityMatrix = eye(nx);
			% Jacobian at of h
			sqdhdx = IdentityMatrix(meas_.states,:);
		elseif (isfield(meas_,'statefcn'))
			% if user provides a measurement state function
			measpred(:,i,k) = meas_.statefcn(x(:,i,k));
			resid(:,i,k) = meas_.data(:,i,k) - meas_.statefcn(x(:,i,k));
			% Jacobian of h at time point i if statefcn provided by user
			sqdhdx = meas_.dstatedx(x(:,i,k));
		elseif (isfield(meas_,'measurematrix'))
			% if output matrix y = measurematrix*x is provided
			measpred(:,i,k) = meas_.measurematrix*(x(:,i,k));
			resid(:,i,k) = meas_.data(:,i,k) - meas_.measurematrix*(x(:,i,k));
			% Jacobian of h at time point i if statefcn provided by user
			sqdhdx = meas_.measurematrix;
		end	

		% squeeze matrix / vector out of struct
		sqsx = sx(:,:,i,k);
		r = resid(:,i,k);
		%% look for bad data (inf, NaN,...) and do not use it for calculation
		badData = find(~isfinite(r));
		r(badData) = 0;
		sqsx(badData,:) = 0;
		% tmp = dh/dp via chain rule
		tmp = sqdhdx*sqsx;
		% update for current time point i and data set k
		phi = phi + r'*meas_.weight*r;
		grad = grad + tmp'*meas_.weight*r;
		Hgn = Hgn + tmp'*meas_.weight*tmp;
	end %i = 1: ntimes
  end %k = 1:nics

  %% gradient with respect to untransformed parameters
	grad = -2*grad';
	obj_.grad = grad;
  %% gradient with respect to log transformed parameters
  %% this is the optimizer's gradient (see gradi function)
	scale = ones(npar,1);
	scale(obj_.logtr) = data.theta(trind);
	obj_.gradlogtr = obj_.grad*diag(scale);
  %% Gauss Newtown approximation of Hessian with untransformed parameters
	obj_.Hgn =  2*Hgn;
	obj_.Hgnlogtr =  2*diag(scale)*Hgn*diag(scale);
	obj_.x = x;
	obj_.sx = sx;
	obj_.resid = resid;
	obj_.measpred = measpred;
  
  % print on screen if printlevel is set accordingly
  if(obj_.printlevel == 1)
	disp('objective function:');
	disp(phi);
	disp('gradient:');
	disp(grad);
  end

%% check on where the optimizer is looking
%theta
%phi
end %Likelihood
%% END OF FUNCTION LIKELIHOOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, sx, status] = cvode_sens(f, fp, x0, estflag, data, sx0, tmeas)
%% Call the cvode integrator and calculate states, x, and 
%% sensitivities, sx, at the requested times.
%% f is rhs of ode, fp is rhs of sensitivity ode
%% status = 0 -> no errors
global mod_

  t0 = 0;
  nics = size(x0,2);

  %% output size
  nx = size(x0,1);
  nt = size(tmeas,1);
  np = size(sx0,2);
  %% output initialization
  x = zeros (nx,nt,nics);
  sx = zeros (nx, np,nt,nics);

  for k = 1:nics %for each IC =  dataset
	
	%% Set options for integrator.
	options = CVodeSetOptions ('UserData', data, ...
				   'RelTol', mod_.rtol, ...
				   'AbsTol', mod_.atol, ...
				   'MaxNumSteps',mod_.odesteps);
	
	%% Allocate storage and set initial conditions for integrator.
%-	CVodeMalloc (f, t0, x0(:,k), options, data);
        CVodeInit(f,'BDF','Newton',t0,x0(:,k),options);
	%% Number of paramters
	np = numel(estflag);

	fsa_options = CVodeSensSetOptions('method','Simultaneous', ...
					 'ErrControl', true,...
					 'ParamField', 'theta',...
					 'ParamList', estflag,...
					 'ParamScales', ones(1, np));

	%% Set options for forward sensitivity problem.
	if (isempty(fp))
	   CVodeSensInit (np, [], sx0(:,:,k), fsa_options);		      
	else
	   CVodeSensInit (np, fp, sx0(:,:,k), fsa_options);		      
	end %if (isempty(fp))

	%% Allocate storage for forward sensitivity problem.
%-	CVodeSensMalloc (np, 'Simultaneous', sx0(:,:,k), fsa_options);

	if (tmeas(1) == t0)  
		% initial condition for this data set k
		x(:,1,k) = x0(:,k);
		sx(:,:,1,k) = sx0(:,:,k);
		iout = 1; 
	else
		iout = 0;
	end

	tout = tmeas(1+iout);
	nsteps = nt;

	while (true)
		[status, t, x_step, sx_step] = CVode (tout, 'Normal');
		if (status == 0)
			iout = iout + 1;
			x(:,iout,k) = x_step;
			sx(:,:,iout,k) = sx_step;
		elseif (status < 0)
			warning ('parest: CVode failed with status = %d', status);
			break;
		end %if (status == 0)
		if (iout == nsteps)
			break;
		end
		tout = tmeas(1 + iout);
	end % while(true)
	%% Free integrator storage.
	CVodeFree ();

  end % for k = 1:nics

end %cvode_sens
%% END OF FUNCTION cvode_sens %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phi = likelihoodA(theta)
%% Calculates negative logarithm of likelihood function,
%% this function is minimized by parest.m.
%% Gradients are calculated by adjoint method.
  global mod_ meas_ obj_ data

  if(obj_.printlevel == 1)
	disp('----------------------');	
	disp('examined parameter values:');
	disp(theta);
  end

  %% data.theta contains all parameters, not only those which are estimated!
  data.theta = mod_.param;
  %% set the model parameters; note the log transformation
  data.theta(obj_.estflag) = theta;
  trind =  obj_.estflag(obj_.logtr);
  data.theta(trind) = exp(theta(obj_.logtr));

  %% set the state ICs that are parameters
  for i = 1:numel(mod_.stateic_par)
	if (mod_.stateic_par(i) > 0)
		mod_.ic(i) =  data.theta(mod_.stateic_par(i));
	end
  end

  t0 = 0;
  nmeas = size(meas_.data, 1);
  ntimes = size(meas_.data, 2);
  nx = size(mod_.ic,1);
  npar = numel(obj_.estflag);

  %% initialize objective and gradient
  phi = 0;
  grad = zeros (1, npar);

  %% Set options for integrator.
  	options = CVodeSetOptions ('UserData', data, ...
				   'RelTol', mod_.rtol, ...
				   'AbsTol', mod_.atol, ...
				   'MaxNumSteps', mod_.odesteps);

  %% number of IC's = data sets
  nics = size(mod_.ic, 2);

  x = zeros (nx,ntimes,nics);
  resid = zeros (nx,ntimes,nics);
  measpred = zeros (nx,ntimes,nics);
  mu = zeros (nx,ntimes,nics);

for k = 1:nics

  %% Allocate storage and set initial conditions for integrator.
  if (~isfield(mod_,'cvodesfcn'))
  CVodeInit (@oderhs, 'BDF', 'Newton', 0, mod_.ic(:,k), options);
  else
  CVodeInit (mod_.cvodesfcn, 'BDF', 'Newton', 0, mod_.ic(:,k), options);
  end

  %% Allocate memory for checkpointing of forward solution for
  %% backward integration of adjoint problem.
  CVodeAdjInit (150, 'Hermite');

  tmeas = meas_.time(:,k);
  t0 = 0;
  t1 = tmeas(end);

  x(:,:,k) = zeros (nx,ntimes);
  mu(:,:,k) = zeros (nx,ntimes);

  %% Integrate forward, saving intermediate points for calculation of
  %% the objective function.
	if (tmeas(1) == t0)  
		% initial condition for this data set k
		x(:,1,k) = mod_.ic(:,k);
		sx(:,:,1,k) = sx0(:,:,k);
		iout = 1; 
	else
		iout = 0;
	end

	tout = tmeas(1+iout);
	nsteps = nt;

	while (true)
		[status, t, x_step, sx_step] = CVode (tout, 'Normal');
		if (status == 0)
			iout = iout + 1;
			x(:,iout,k) = x_step;
		elseif (status < 0)
			warning ('parest: CVode failed with status = %d', status);
			CVodeFree ();
			phi = realmax;
			return;
		end %if (status == 0)
		if (iout == nsteps)
			break;
		end
		tout = tmeas(1 + iout);
	end % while(true)

  if(status < 0)
	CVodeFree ();
	phi = realmax;
	return;
  end

  %%compute residual and objective function for these param values
  %%compute gradient, here NO Hessian approximation; have to store in common
  if (isfield(meas_,'states'))
		% if user provides a vector of measured states
		statefcnExists = false;
		measpred(:,:,k) = x(meas_.states,:,k);
		resid(:,:,k) = meas_.data(:,:,k) - x(meas_.states,:,k);
		IdentityMatrix = eye(nx);
		dhdx = IdentityMatrix(meas_.states,:);
	elseif (isfield(meas_,'measurematrix'))
		% if user provides measurement matrix y = measurematrix*x
		statefcnExists = false;
		measpred(:,:,k) = meas_.measurematrix*x(:,:,k);
		resid(:,:,k) = meas_.data(:,:,k) - meas_.measurematrix*x(:,:,k);
		dhdx = meas_.measurematrix;
	else
		% if user provides a measurement state function
		statefcnExists = true;
		for t_ind = 1:ntimes	%loop over all time points
			resid(:,t_ind,k) = meas_.data(:,t_ind,k) - meas_.statefcn(x(:,t_ind,k));
			measpred(:,t_ind,k) = meas_.statefcn(x(:,t_ind,k));
		end
  end

  obj_.resid = resid;
  obj_.measpred = measpred;
  %% find and delete bad data
  resid(find(~isfinite(resid))) = 0;
	
  %% update objective function
  phi = phi + trace(resid(:,:,k)'*meas_.weight*resid(:,:,k));

  %% Free storage for the integration.
  %CVodeFree ();

  %% Set up reverse integration.  Must be done after the forward
  %% integration is complete.

  %% Initial condition for the backward quadrature.
  w = zeros (1, npar);

  %% Set options for the backward integration of the adjoint equation and
  %% quadrature.
  optionsb = CVodeSetOptions ('RelTol', mod_.rtol, ...
			      'AbsTol', mod_.atol);

  optionsQB = CVodeQuadSetOptions('ErrControl', true, ...
				  'RelTol', sqrt(eps), ...
				  'AbsTol', sqrt(eps));



  %% Initial condition for the backward integration of the adjoint
  %% equations
  if (statefcnExists)
	% Jacobian at last time point if statefcn provided by user
	sqdhdx = meas_.dstatedx(x(:,end,k));
	%% find and delete bad data
	sqdhdx(find(~isfinite(meas_.data(:,end,k))),:) = 0;
  else
	% Jacobian at last time point if only measurement provided by user
	sqdhdx = dhdx;
	%% find and delete bad data
	sqdhdx(find(~isfinite(meas_.data(:,end,k))),:) = 0;
  end
  
  mu_init = -2*( resid(:,end,k)'*meas_.weight*sqdhdx )';
  mu(:,end,k) = mu_init;
  %% Allocate storage and set initial conditions for integrator.
  idxB = CVodeInitB (@adjrhs, 'BDF', 'Newton', t1, mu_init, optionsb);

  CVodeQuadInitB(idxB, @sensrhsA, w, optionsQB);
  %% Integrate backward, accumulating in w the gradient of phi, which is
  %% the least-squares objective \sum resid*meas_.weight*resid.
  for (i = ntimes-1:-1:1)
	tinit = tmeas(i+1);
	tout = tmeas(i);
	CVodeReInitB (idxB, tinit, mu_init, optionsb);

   	while (true)
		[status, t_step, mu_step, w_step] = CVodeB (tout, 'Normal');

		if (status < 0)
			warning ('parest: CVodeB failed with status = %d', status);
			CVodeFree ();
			phi = realmax;
			return;
		end

		%% We are integrating backward, and we are storing the results in
		%% reverse order too.
    		if (status == 0 || (status == 1 && t_step == 0))
		if (statefcnExists)
			% Jacobian at time point i if statefcn provided by user
			sqdhdx = meas_.dstatedx(x(:,i,k));
			%% find and delete bad data
			sqdhdx(find(~isfinite(meas_.data(:,i,k))),:) = 0;
		else
			% Jacobian at time point i if no measurement fcn provided by user
			sqdhdx = dhdx;
			%% find and delete bad data
			sqdhdx(find(~isfinite(meas_.data(:,i,k))),:) = 0;
		end
			mu_init = mu_step - 2*(resid(:,i,k)'*meas_.weight*sqdhdx )';
			mu(:,i,k) = mu_init;
			w = w + w_step';
			break;
		end % if (status == ...
	end % while(true)
  end % for i = ntimes:...

  if( tmeas(1) > t0) %% integrate till t = 0
	CVodeReInitB (idxB, tmeas(1), mu_init, optionsb);
	[status, t_s, mu_0, w_0] = CVodeB (t0, 'Normal');
	if (status ~= 0)
		warning ('parest: CVodeB failed with status = %d', status);
		return;
	else
		w = w + w_step';
   	end
  else
  	mu_0 = mu_init;
  end

  %% Free storage for the integration.
	CVodeFree ();
  %CVadjFree ();

    %% gradient with respect to untransformed parameters
    %% Eqn (3.19) of CVodeS user guide
    %% sx0 = 0 for parameters that are not state IC
    %% sx0 = 1 for parameters that are state IC
    grad = grad + w + mu_0'*mod_.sx0(:,:,k);

end % for k = 1:nics

    obj_.grad = grad;
    %% gradient with respect to log transformed parameters
    %% this is the optimizer's gradient (see gradi function)
    scale = ones(npar,1);
    scale(trind) = data.theta(trind);
    obj_.gradlogtr = obj_.grad*diag(scale);
    %% values of x at tmeas
    obj_.x = x;
    obj_.mu = mu;

  % print on screen if printlevel is set accordingly
  if(obj_.printlevel == 1)
	disp('objective function:');
	disp(phi);
	disp('gradient:');
	disp(grad);
  end

end %LikelihoodA
%% END OF FUNCTION LIKELIHOODA.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = hineq(theta)
%% provides h_ineq(x) for inequality constraint h_ineq(x) >= 0
%% to ensure that parameters are between lower and upper bound
%% parlb <= theta <= parub
  global obj_
  retval = [theta - obj_.parlb; -theta + obj_.parub];
end %hineq
%% END OF FUNCTION hineq %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function retval = gradi(theta)
%% optimizer should have just called the function likelihood with the same
%% theta value, and likelihood has computed obj_.gradlogtr.
%% So pull the gradient out of the global variables and return
  global obj_ likelihood_last_theta mod_ solveroption

	if (~ isequal (theta, likelihood_last_theta))
		if(solveroption == 3)
		%% resolve model to update cache
		obj = likelihoodA (theta);
		else
		obj = likelihood (theta);
		end
	end
  retval = obj_.gradlogtr';
end %gradi
%% END OF FUNCTION gradi %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xdot, flag, new_data] = oderhs(t, x, data)
%% Right hand side of differential equation
%% dx/dt = f(x,t,theta)
  global mod_

  [xdot] = mod_.odefcn(x, t, data.theta);
  flag = 0;
  new_data = [];
end %oderhs
%% END OF FUNCTION oderhs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xsd, flag, new_data] = sensrhs(t, x, xdot, xs, data)
%% Right hand side of forward sensitivity equation
%% sx = dx/dp
%% dsx / dt = df/dx * sx + df/dp
  global mod_ obj_

  dfdp = mod_.dodedp(x, t, data.theta);
  dfdp = dfdp(:,obj_.estflag);
  xsd = mod_.dodedx(x, t, data.theta)*xs + dfdp;
  flag = 0;
  new_data = [];
end %sensrhs
%% END OF FUNCTION sensrhs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qbd, flag, new_data] = sensrhsA (t, x, mu, data)
%% Backward quadrature function use to compute gradient of objective
%% function.  In this problem, it computes Equation 3.19 from the CVodes
%% user guide:
%%   dg/dp(t1) = mu(t0)s(t0) + g_p(t1) + \int_t0^t1 mu'*f_p
%% 
%% compute the integral by integrating
%%    dw/dt = - mu' f_p
%%    w(t1) = 0
%% 
%%  backward over the time horizon.
  global mod_ obj_

  dfdp = mod_.dodedp(x, t, data.theta);
  dfdp = dfdp(:,obj_.estflag);
  qbd = -(mu'*dfdp)';
  flag = 0;
  new_data = [];
end %sensrhsA
%% END OF FUNCTION sensrhsA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mudot, flag, new_data] = adjrhs (t, x, mu, data)
%% RHS of adjoint system.  In this problem, it computes Equation 3.20
%% from the CVodes user guide:
%%    mudot = - f_x'*mu
%%    mu(t1) = g_x' at t=t1
  global mod_

  mudot = -mod_.dodedx(x, t, data.theta)'*mu;
  flag = 0;
  new_data = [];
end %adjrhs
%% END OF FUNCTION adjrhs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

