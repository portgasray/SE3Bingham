function expression = cBinghamGradLogNorm6(in1)
%CBINGHAMGRADLOGNORM6
%    EXPRESSION = CBINGHAMGRADLOGNORM6(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    13-Mar-2019 17:55:25

x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
x5 = in1(5,:);
x6 = in1(6,:);
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
t22 = x1-x6;
t23 = 1.0./t22;
t24 = x2-x6;
t25 = 1.0./t24;
t26 = x3-x6;
t27 = 1.0./t26;
t28 = x4-x6;
t29 = 1.0./t28;
t30 = x5-x6;
t31 = 1.0./t30;
t32 = exp(x1);
t33 = t3.*t5.*t9.*t15.*t23.*t32;
t34 = exp(x2);
t35 = 1.0./t2.^2;
t36 = exp(x3);
t37 = 1.0./t4.^2;
t38 = exp(x4);
t39 = 1.0./t8.^2;
t40 = exp(x5);
t41 = 1.0./t14.^2;
t42 = exp(x6);
t43 = 1.0./t22.^2;
t44 = t5.*t7.*t13.*t19.*t27.*t36;
t45 = t15.*t17.*t19.*t21.*t31.*t40;
t49 = t3.*t7.*t11.*t17.*t25.*t34;
t54 = t9.*t11.*t13.*t21.*t29.*t38;
t55 = t23.*t25.*t27.*t29.*t31.*t42;
t46 = t33+t44+t45-t49-t54-t55;
t47 = 1.0./t46;
t48 = t5.*t9.*t15.*t23.*t32.*t35;
t50 = 1.0./t6.^2;
t51 = 1.0./t10.^2;
t52 = 1.0./t16.^2;
t53 = 1.0./t24.^2;
t56 = t3.*t9.*t15.*t23.*t32.*t37;
t57 = t3.*t11.*t17.*t25.*t34.*t50;
t58 = t7.*t13.*t19.*t27.*t36.*t37;
t59 = 1.0./t12.^2;
t60 = 1.0./t18.^2;
t61 = 1.0./t26.^2;
t62 = t3.*t5.*t15.*t23.*t32.*t39;
t63 = t3.*t7.*t17.*t25.*t34.*t51;
t64 = t9.*t11.*t21.*t29.*t38.*t59;
t65 = t9.*t13.*t21.*t29.*t38.*t51;
t66 = 1.0./t20.^2;
t67 = 1.0./t28.^2;
t68 = t3.*t5.*t9.*t23.*t32.*t41;
t69 = t3.*t7.*t11.*t25.*t34.*t52;
t70 = t15.*t17.*t19.*t31.*t40.*t66;
t71 = t17.*t19.*t21.*t31.*t40.*t41;
t72 = 1.0./t30.^2;
t73 = t3.*t5.*t9.*t15.*t32.*t43;
t74 = t3.*t7.*t11.*t17.*t34.*t53;
t75 = t23.*t25.*t27.*t29.*t42.*t72;
t76 = t23.*t25.*t29.*t31.*t42.*t61;
t77 = t23.*t27.*t29.*t31.*t42.*t53;
expression = [-t47.*(-t33+t48+t56+t58+t62+t68+t71+t73-t7.*t11.*t17.*t25.*t34.*t35-t11.*t13.*t21.*t29.*t38.*t39-t25.*t27.*t29.*t31.*t42.*t43),t47.*(t48-t49+t57+t63+t65+t69+t74+t77-t7.*t11.*t17.*t25.*t34.*t35-t5.*t13.*t19.*t27.*t36.*t50-t15.*t19.*t21.*t31.*t40.*t52),t47.*(t44+t56-t57+t58+t64+t76-t5.*t7.*t13.*t19.*t36.*t61-t5.*t7.*t13.*t27.*t36.*t60+t5.*t13.*t19.*t27.*t36.*t50-t5.*t7.*t19.*t27.*t36.*t59-t15.*t17.*t21.*t31.*t40.*t60),-t47.*(t54-t62+t63+t64+t65+t70+t11.*t13.*t21.*t29.*t38.*t39-t5.*t7.*t19.*t27.*t36.*t59-t9.*t11.*t13.*t21.*t38.*t67-t9.*t11.*t13.*t29.*t38.*t66-t23.*t25.*t27.*t31.*t42.*t67),t47.*(t45+t68-t69+t70+t71+t75+t5.*t7.*t13.*t27.*t36.*t60-t9.*t11.*t13.*t29.*t38.*t66+t15.*t19.*t21.*t31.*t40.*t52-t15.*t17.*t19.*t21.*t40.*t72+t15.*t17.*t21.*t31.*t40.*t60),-t47.*(t55-t73+t74+t75+t76+t77-t5.*t7.*t13.*t19.*t36.*t61+t9.*t11.*t13.*t21.*t38.*t67-t15.*t17.*t19.*t21.*t40.*t72+t25.*t27.*t29.*t31.*t42.*t43+t23.*t25.*t27.*t31.*t42.*t67)];
