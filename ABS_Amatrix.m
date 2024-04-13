function A = ABS_Amatrix(Cd,J,R,m,omega,slope,theta1,theta2,theta3,v,wind)
%ABS_Amatrix
%    A = ABS_Amatrix(Cd,J,R,M,OMEGA,SLOPE,THETA1,THETA2,THETA3,V,WIND)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    30-Apr-2023 16:38:00

t2 = cos(slope);
t3 = abs(v);
t4 = sign(v);
t5 = R.*omega;
t6 = 1.0./J;
t7 = 1.0./v;
t8 = t7.^2;
t9 = -t5;
t10 = 1.0./t3;
t12 = t7.*theta3;
t11 = t10.^2;
t13 = t9+v;
t14 = R.*t12;
t17 = -t12.*(t5-v);
t18 = -t8.*theta3.*(t5-v);
t21 = t8.*theta3.*(t5-v);
t15 = abs(t13);
t16 = sign(t13);
t19 = t10.*t15.*theta2;
t20 = t10.*t16.*theta2;
t22 = t4.*t11.*t15.*theta2;
t23 = -t19;
t25 = -t22;
t24 = exp(t23);
t28 = t20+t25;
t26 = -t24;
t27 = R.*t20.*t24;
t31 = t24.*t28;
t29 = t17+t26+1.0;
t30 = t14+t27;
t32 = t12+t21+t31;
mt1 = [-(Cd.*(v.*2.0-wind.*2.0))./m+t2.*t32.*theta1.*(9.81e+2./1.0e+2),R.*m.*t2.*t6.*t32.*theta1.*(-9.81e+2./1.0e+2),0.0,0.0,0.0,t2.*t30.*theta1.*(-9.81e+2./1.0e+2),R.*m.*t2.*t6.*t30.*theta1.*(9.81e+2./1.0e+2),0.0,0.0,0.0,t2.*(t24+t12.*(t5-v)-1.0).*(-9.81e+2./1.0e+2),R.*m.*t2.*t6.*(t24+t12.*(t5-v)-1.0).*(9.81e+2./1.0e+2),0.0,0.0,0.0,t2.*t10.*t15.*t24.*theta1.*(9.81e+2./1.0e+2),R.*m.*t2.*t6.*t10.*t15.*t24.*theta1.*(-9.81e+2./1.0e+2),0.0,0.0,0.0,t2.*t7.*theta1.*(t5-v).*(-9.81e+2./1.0e+2)];
mt2 = [R.*m.*t2.*t6.*t7.*theta1.*(t5-v).*(9.81e+2./1.0e+2),0.0,0.0,0.0];
A = reshape([mt1,mt2],5,5);
