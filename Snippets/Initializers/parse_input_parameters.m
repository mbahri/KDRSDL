function [ params ] = parse_input_parameters( arguments, params )
%PARSE_INPUT_PARAMETERS Parse the input params. and sets the structure
%
% Mehdi Bahri - Imperial College London
% July, 2016
%
% Last modified August, 2017

if (rem(length(arguments),2)==1)
  error('Optional parameters should always go by pairs');
else
  for i=1:2:(length(arguments)-1)
    if ~ischar (arguments{i})
      error (['Unknown type of optional parameter name (parameter' ...
	      ' names must be strings).']);
    end
    % change the value of parameter
    switch lower (arguments{i})
     case 'r'
      params.r = arguments{i+1};
     case 'lambda'
      params.lambda = arguments{i+1};
     case 'tol'
      params.tol = arguments{i+1};
     case 'maxiter'
      params.MAXITER = arguments{i+1};
     case 'rho'
      params.rho = arguments{i+1};
     case 'enforce_monotonicity'
      params.MONOTONE = arguments{i+1};
     case 'alpha_a'
      params.alpha_a = arguments{i+1};
     case 'alpha_b'
      params.alpha_b = arguments{i+1};
     case 'alpha'
      params.alpha = arguments{i+1};
     case 'mu'
      params.mu = arguments{i+1};
     case 'ground_o'
      params.g_O = arguments{i+1};
     case 'ground_e'
      params.g_E = arguments{i+1};
     case 'mean'
      params.MEAN = arguments{i+1};
     case 'parallel'
      params.PARALLEL = arguments{i+1};
     case 'time'
      params.TIME = arguments{i+1};
     otherwise
      error(['Unrecognized parameter: ''' arguments{i} '''']);
    end
  end
end

end

