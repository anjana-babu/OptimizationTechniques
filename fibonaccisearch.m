function [val,xmin]=fibonaccisearch(a,b,min_interval,function_name,max_iter)
% FIBONACCISEARCH A linear search procedure for minimizing a strictly
% quasiconvex function over a closed bounded interval.
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
%   [val,x] = fibonaccisearch(0,2,0.02,'0.5*x^2-sin(x)')
%   [val,x] = fibonaccisearch(0,2,0.02,'0.5*x^2-sin(x)',25)
%   [val,x] = fibonaccisearch(-3,5,0.002,'x^2+2*x',50)
%
% Developed by Anjana Babu on 27-Apr-2012

if (nargin < 5)
    max_iter = 100; % Setting maximum iterations to 100 if not specified by user
end

f = inline(function_name); % Getting the function

n = calculatesteps(a,b,min_interval); % Calculating number of steps n needed

% Initialization Step
lambda = a + (fibnum(n-2)/fibnum(n))*(b - a);
mu = a + (fibnum(n-1)/fibnum(n))*(b - a);

f_lambda = f(lambda);
f_mu = f(mu);

% Plotting
figure; hold on;
plot(lambda,f_lambda,'rx');
plot(mu,f_mu,'rx');


% Main Step
k=1;

while (((b-a) > min_interval) && (k < max_iter))   
    if(f_lambda>f_mu) % Interval [lambda,b] is chosen
        a = lambda;
        lambda = mu;
        mu = a + (fibnum(n-k-1)/fibnum(n-k))*(b - a);
        f_lambda = f_mu;
        f_mu = f(mu);
        plot(mu,f_mu,'rx');
    else  % Interval [a,mu] is chosen
        b = mu;
        mu = lambda;
        lambda = a + (fibnum(n-k-2)/fibnum(n-k))*(b - a);
        f_mu = f_lambda;
        f_lambda = f(lambda);
        plot(lambda,f_lambda,'rx');
    end % End of if construct
    k = k+1; % Incrementing k
end % End of while loop

if(f_lambda>f_mu)
    a = lambda;
else
    b = mu;
end

xmin = (a+b)/2; % Setting output variable xmin
val=f(xmin); % Setting output variable out
plot(xmin,val,'b*');

end % End of function fibonaccisearch




function n = calculatesteps(a_init,b_init,l)
% CALCULATESTEPS Calculates no of steps needed in the fibonacci search
Fn = (b_init-a_init)/l;
i=0;
while (fibnum(i)<Fn)
    i = i+1;
end
n = i;
end




function f = fibnum(n)
% FIBNUM Fibonacci number.
% FIBNUM(n) generates the nth Fibonacci number.
if n <= 2
f = 1;
else
f = fibnum(n-1) + fibnum(n-2);
end
end % End of function fibnum

