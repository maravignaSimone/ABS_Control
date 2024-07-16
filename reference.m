function r = reference(t, theta1, theta2, theta3)
    
mu = @(lambda) sign(lambda) * theta1 * (1-exp(-abs(lambda)*theta2))-(lambda*theta3);
[lambda_min, ~] = fminsearch(mu, 0);

r = lambda_min;
end

