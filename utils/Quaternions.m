clear all;
syms d0 d1 d2 d3 d4 d5 d6 d7 real
r = [d0 d1 d2 d3];
d = [d4 d5 d6 d7];

r0 = d0
d0 = d4
vec_r = r(2:4);
vec_d = d(2:4);

inv_vec_d = -vec_d;

% (d, r*) multiplication of two quaternions
quat_pos_sca = d0*r0 - dot(inv_vec_d, vec_r);
quat_pos_vec = d0.*vec_r + r0.*inv_vec_d + cross(inv_vec_d, vec_r);
% not quat_pos_sca + quat_pos_vec
quat_pos = [quat_pos_sca quat_pos_vec];

% sanity check
conju_r = [d0, -d1, -d2, -d3];
quatmultiply(d, conju_r) - quat_pos

% for Hamilton operator is invertible?
% syms dx dy dz real
% Hs = [0, -dx, -dy, -dz; dx, 0, -dz, dy; dy, dz, 0, -dx; dz, -dy, dx, 0];
% inv(Hs) %if (dx^2 + dy^2 + dz^2) != 0