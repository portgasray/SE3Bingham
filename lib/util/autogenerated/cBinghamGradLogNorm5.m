function expression = cBinghamGradLogNorm5(in1)
%CBINGHAMGRADLOGNORM5
%    EXPRESSION = CBINGHAMGRADLOGNORM5(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    13-Mar-2019 17:55:22

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
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
t14 = x1-x5;
t15 = 1.0./t14;
t16 = x2-x5;
t17 = 1.0./t16;
t18 = x3-x5;
t19 = 1.0./t18;
t20 = x4-x5;
t21 = 1.0./t20;
t22 = exp(x1);
t23 = t3.*t5.*t9.*t15.*t22;
t24 = exp(x2);
t25 = 1.0./t2.^2;
t26 = exp(x3);
t27 = 1.0./t4.^2;
t28 = exp(x4);
t29 = 1.0./t8.^2;
t30 = exp(x5);
t31 = 1.0./t14.^2;
t32 = t5.*t7.*t13.*t19.*t26;
t33 = t15.*t17.*t19.*t21.*t30;
t37 = t3.*t7.*t11.*t17.*t24;
t41 = t9.*t11.*t13.*t21.*t28;
t34 = t23+t32+t33-t37-t41;
t35 = 1.0./t34;
t36 = t5.*t9.*t15.*t22.*t25;
t38 = 1.0./t6.^2;
t39 = 1.0./t10.^2;
t40 = 1.0./t16.^2;
t42 = t3.*t9.*t15.*t22.*t27;
t43 = t3.*t11.*t17.*t24.*t38;
t44 = t7.*t13.*t19.*t26.*t27;
t45 = 1.0./t12.^2;
t46 = 1.0./t18.^2;
t47 = t3.*t5.*t15.*t22.*t29;
t48 = t3.*t7.*t17.*t24.*t39;
t49 = t9.*t11.*t21.*t28.*t45;
t50 = t9.*t13.*t21.*t28.*t39;
t51 = 1.0./t20.^2;
t52 = t3.*t5.*t9.*t22.*t31;
t53 = t3.*t7.*t11.*t24.*t40;
t54 = t15.*t17.*t19.*t30.*t51;
t55 = t17.*t19.*t21.*t30.*t31;
expression = [-t35.*(-t23+t36+t42+t44+t47+t52+t55-t7.*t11.*t17.*t24.*t25-t11.*t13.*t21.*t28.*t29),t35.*(t36-t37+t43+t48+t50+t53-t7.*t11.*t17.*t24.*t25-t5.*t13.*t19.*t26.*t38-t15.*t19.*t21.*t30.*t40),t35.*(t32+t42-t43+t44+t49-t5.*t7.*t13.*t26.*t46+t5.*t13.*t19.*t26.*t38-t5.*t7.*t19.*t26.*t45-t15.*t17.*t21.*t30.*t46),-t35.*(t41-t47+t48+t49+t50+t54-t5.*t7.*t19.*t26.*t45+t11.*t13.*t21.*t28.*t29-t9.*t11.*t13.*t28.*t51),t35.*(t33+t52-t53+t54+t55+t5.*t7.*t13.*t26.*t46-t9.*t11.*t13.*t28.*t51+t15.*t19.*t21.*t30.*t40+t15.*t17.*t21.*t30.*t46)];
