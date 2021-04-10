clear all
syms v0 v1 v2 v3 w0 w1 w2 w3 real
v = [v0 v1 v2 v3];
w = [w0 w1 w2 w3];

v_ = v(2:4); %broke in vector components
w_ = w(2:4);
r_scalar = v0*w0-dot(v_,w_)

r_vector = v0.*w_ + v_*w0 + cross(v_,w_)

r = [r_scalar r_vector]
quatmultiply(v,w)-r