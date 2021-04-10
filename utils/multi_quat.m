clear all;
syms d0 d1 d2 d3 d4 d5 d6 d7 real
v = [d0 d1 d2 d3];
w = [d4 d5 d6 d7];

conjugate_v = [d0 -d1 -d2 -d3];
conjugate_w = [d4 -d5 -d6 -d7];
v0 = d0;
w0 = d4;
v_ = -v(2:4);
w_ =  w(2:4);

translation_scalar = v0*w0-dot(v_,w_)

translation_vector = v0.*w_ + v_*w0 + cross(v_,w_)

translation = [translation_scalar translation_vector]
%quatmultiply(conjugate_v,w)-r
translation

%norm of q_d: (q_d, q_d*)
%quatmultiply