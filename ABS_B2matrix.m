function B2 = ABS_B2matrix(Cd,J,R,m,omega,slope,theta1,theta2,theta3,v,wind)
%ABS_B2matrix
%    B2 = ABS_B2matrix(Cd,J,R,M,OMEGA,SLOPE,THETA1,THETA2,THETA3,V,WIND)

%    This function was generated by the Symbolic Math Toolbox version 23.2.
%    25-May-2024 10:30:34

t2 = abs(v);
t3 = sin(slope);
t4 = R.*omega;
t5 = 1.0./v;
t6 = -t4;
t7 = 1.0./t2;
t10 = -t5.*(t4-v);
t8 = t6+v;
t11 = sign(t10);
t12 = t10.*theta3;
t9 = abs(t8);
t13 = t7.*t9.*theta2;
t14 = -t13;
t15 = exp(t14);
t16 = t15-1.0;
t17 = t11.*t16.*theta1;
t18 = t12+t17;
B2 = reshape([(Cd.*(v.*2.0-wind.*2.0))./m,0.0,0.0,0.0,0.0,cos(slope).*(-9.81e+2./1.0e+2)-t3.*t18.*(9.81e+2./1.0e+2),(R.*m.*t3.*t18.*(9.81e+2./1.0e+2))./J,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[5,5]);
end
