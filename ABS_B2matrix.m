function B2 = ABS_B2matrix(Cd,J,R,m,omega,slope,theta1,theta2,theta3,v,wind)
%ABS_B2matrix
%    B2 = ABS_B2matrix(Cd,J,R,M,OMEGA,SLOPE,THETA1,THETA2,THETA3,V,WIND)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    30-Apr-2023 16:38:00

t2 = abs(v);
t3 = sin(slope);
t4 = R.*omega;
t5 = 1.0./v;
t6 = -t4;
t7 = 1.0./t2;
t10 = -t5.*theta3.*(t4-v);
t8 = t6+v;
t9 = abs(t8);
t11 = t7.*t9.*theta2;
t12 = -t11;
t13 = exp(t12);
t14 = -t13;
B2 = reshape([(Cd.*(v.*2.0-wind.*2.0))./m,0.0,0.0,0.0,0.0,cos(slope).*(9.81e+2./1.0e+2)+t3.*theta1.*(t13+t5.*theta3.*(t4-v)-1.0).*(9.81e+2./1.0e+2),(R.*m.*t3.*theta1.*(t13+t5.*theta3.*(t4-v)-1.0).*(-9.81e+2./1.0e+2))./J,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[5,4]);
