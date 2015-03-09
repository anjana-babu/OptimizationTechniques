function [val,xmin]=goldensection(a,b,min_interval,function_name,max_iter)
% GOLDENSECTION to find the minimum of a unimodal continuous function over an interval without using derivatives
%   Inputs:
%   ------
%   'a'                 is the starting of the initial interval
%   'b'                 is the ending of the initial interval
%   'min_interval'      is the maximum length allowed for the final
%   interval of uncertainty
%   'function_name'     is the function for which the min value has to be
%   found. Eg, '0.5*x^2-2*sin(x)'
%   'max_iter'          is the maximum no of iterations that can be taken 
%   Outputs:
%   -------
%   'val'               is the minimum value of the function in the given
%   interval
%   'xmin'              is the x at which function becomes minimum
%   Examples:
%   --------
%   [val,x] = goldensection(0,2,0.02,'0.5*x^2-sin(x)')
%   [val,x] = goldensection(0,2,0.02,'0.5*x^2-sin(x)',25)
%
% Developed by Anjana Babu on 27-Apr-2012

if (nargin < 5)
    max_iter = 100; % Setting maximum iterations to 100 if not specified by user
end

f = inline(function_name); % Getting the function

alpha = 0.618;

% Initialization Step
lambda = a + (1 - alpha)*(b - a);
mu = a + (alpha)*(b - a);

f_lambda = f(lambda);
f_mu = f(mu);

% Plotting
figure; hold on;
plot(lambda,f_lambda,'rx');
plot(mu,f_mu,'rx');

% Main Step
k=0;
while (((b-a) > min_interval) && (k < max_iter))
    if(f_lambda>f_mu) % Interval [lambda,b] is chosen
        a = lambda;
        lambda = mu;
        mu = a + (alpha)*(b - a);
        f_lambda = f_mu;
        f_mu = f(mu);
        plot(mu,f_mu,'rx');
    else  % Interval [a,mu] is chosen
        b = mu;
        mu = lambda;
        lambda = a + (1 - alpha)*(b - a);
        f_mu = f_lambda;
        f_lambda = f(lambda);
        plot(lambda,f_lambda,'rx');
    end % End of if construct
    
    k = k+1; % Incrementing k
end % End of while loop

xmin = (a+b)/2; % Setting output variable xmin
val=f(xmin); % Setting output variable out
plot(xmin,val,'b*');

end % End of function goldensection


