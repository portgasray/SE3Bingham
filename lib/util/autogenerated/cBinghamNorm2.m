function expression = cBinghamNorm2(in1)
%CBINGHAMNORM2
%    EXPRESSION = CBINGHAMNORM2(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    13-Mar-2019 17:55:18

x1 = in1(1,:);
x2 = in1(2,:);
t2 = x1-x2;
t3 = 1.0./t2;
expression = pi.^2.*(t3.*exp(x1)-t3.*exp(x2)).*2.0;