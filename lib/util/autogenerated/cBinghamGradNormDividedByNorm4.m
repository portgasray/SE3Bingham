function expression = cBinghamGradNormDividedByNorm4(in1)
%CBINGHAMGRADNORMDIVIDEDBYNORM4
%    EXPRESSION = CBINGHAMGRADNORMDIVIDEDBYNORM4(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    13-Mar-2019 17:55:21

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
t2 = x1-x2;
t3 = 1.0./t2;
t4 = x1-x3;
t5 = 1.0./t4;
t6 = x2-x3;
t7 = 1.0./t6;
t8 = x1-x4;
t9 = 1.0./t8;
t10 = x2-x4;
t11 = 1.0./t10;
t12 = x3-x4;
t13 = 1.0./t12;
t14 = exp(x1);
t15 = t3.*t5.*t9.*t14;
t16 = exp(x2);
t17 = 1.0./t2.^2;
t18 = exp(x3);
t19 = 1.0./t4.^2;
t20 = exp(x4);
t21 = 1.0./t8.^2;
t22 = t5.*t7.*t13.*t18;
t26 = t3.*t7.*t11.*t16;
t29 = t9.*t11.*t13.*t20;
t23 = t15+t22-t26-t29;
t24 = 1.0./t23;
t25 = t5.*t9.*t14.*t17;
t27 = 1.0./t6.^2;
t28 = 1.0./t10.^2;
t30 = t3.*t9.*t14.*t19;
t31 = t3.*t11.*t16.*t27;
t32 = t7.*t13.*t18.*t19;
t33 = 1.0./t12.^2;
t34 = t3.*t5.*t14.*t21;
t35 = t3.*t7.*t16.*t28;
t36 = t9.*t11.*t20.*t33;
t37 = t9.*t13.*t20.*t28;
expression = [-t24.*(-t15+t25+t30+t32+t34-t7.*t11.*t16.*t17-t11.*t13.*t20.*t21),t24.*(t25-t26+t31+t35+t37-t7.*t11.*t16.*t17-t5.*t13.*t18.*t27),t24.*(t22+t30-t31+t32+t36-t5.*t7.*t18.*t33+t5.*t13.*t18.*t27),-t24.*(t29-t34+t35+t36+t37-t5.*t7.*t18.*t33+t11.*t13.*t20.*t21)];
