syms a0 a1 a2 a3 a4 a5 l1 l2 l3 l4 l5 l6 rx ry rz px py pz;
syms j00 j01 j02 j03 j04 j05 
syms j10 j11 j12 j13 j14 j15 
syms j20 j21 j22 j23 j24 j25 
syms j30 j31 j32 j33 j34 j35 
syms j40 j41 j42 j43 j44 j45 
syms j50 j51 j52 j53 j54 j55 

t = @(x)([-cos(x) sin(x) 0;0 0 1;sin(x) cos(x) 0]);
r = @(l)([0;0;l]);
td = @(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0; 0 0 1]);

r = t(a0) * r(l1) + t(a0)*t(a1) * r(l2) + t(a0)*t(a1)*t(a2) * r(l3) + t(a0)*t(a1)*t(a2)*t(a3) * r(l4) + t(a0)*t(a1)*t(a2)*t(a3)*t(a4) * r(l5) + t(a0)*t(a1)*t(a2)*t(a3)*t(a4)*t(a5) * r(l6);
p = (td(a0) * td(a1) * td(a2) * td(a3) * td(a4) * td(a5)) * [0;0;1];
f = [r;p];

% #include <stdio.h>
% int main(void)
% {
%     int i, j;
%     for (i=0; i<6; i++)
%     {
%         for (j=0; j < 6; j++)
%         {
%             printf("j%d%d = diff(f(%d,1), a%d);\n", i, j, i+1, j);
%         }
%     }

%     for (i=0; i<6; i++)
%     {
%         for (j=0; j < 6; j++)
%         {
%             printf("fprintf(\"jacobian[%d][%d] = %%s;\\n\", j%d%d);\n", i, j, i, j);
%         }
%     }
%     return 0;
% }
j00 = diff(f(1,1), a0);
j01 = diff(f(1,1), a1);
j02 = diff(f(1,1), a2);
j03 = diff(f(1,1), a3);
j04 = diff(f(1,1), a4);
j05 = diff(f(1,1), a5);
j10 = diff(f(2,1), a0);
j11 = diff(f(2,1), a1);
j12 = diff(f(2,1), a2);
j13 = diff(f(2,1), a3);
j14 = diff(f(2,1), a4);
j15 = diff(f(2,1), a5);
j20 = diff(f(3,1), a0);
j21 = diff(f(3,1), a1);
j22 = diff(f(3,1), a2);
j23 = diff(f(3,1), a3);
j24 = diff(f(3,1), a4);
j25 = diff(f(3,1), a5);
j30 = diff(f(4,1), a0);
j31 = diff(f(4,1), a1);
j32 = diff(f(4,1), a2);
j33 = diff(f(4,1), a3);
j34 = diff(f(4,1), a4);
j35 = diff(f(4,1), a5);
j40 = diff(f(5,1), a0);
j41 = diff(f(5,1), a1);
j42 = diff(f(5,1), a2);
j43 = diff(f(5,1), a3);
j44 = diff(f(5,1), a4);
j45 = diff(f(5,1), a5);
j50 = diff(f(6,1), a0);
j51 = diff(f(6,1), a1);
j52 = diff(f(6,1), a2);
j53 = diff(f(6,1), a3);
j54 = diff(f(6,1), a4);
j55 = diff(f(6,1), a5);
fprintf("jacobian[0][0] = %s;\n", j00);
fprintf("jacobian[0][1] = %s;\n", j01);
fprintf("jacobian[0][2] = %s;\n", j02);
fprintf("jacobian[0][3] = %s;\n", j03);
fprintf("jacobian[0][4] = %s;\n", j04);
fprintf("jacobian[0][5] = %s;\n", j05);
fprintf("jacobian[1][0] = %s;\n", j10);
fprintf("jacobian[1][1] = %s;\n", j11);
fprintf("jacobian[1][2] = %s;\n", j12);
fprintf("jacobian[1][3] = %s;\n", j13);
fprintf("jacobian[1][4] = %s;\n", j14);
fprintf("jacobian[1][5] = %s;\n", j15);
fprintf("jacobian[2][0] = %s;\n", j20);
fprintf("jacobian[2][1] = %s;\n", j21);
fprintf("jacobian[2][2] = %s;\n", j22);
fprintf("jacobian[2][3] = %s;\n", j23);
fprintf("jacobian[2][4] = %s;\n", j24);
fprintf("jacobian[2][5] = %s;\n", j25);
fprintf("jacobian[3][0] = %s;\n", j30);
fprintf("jacobian[3][1] = %s;\n", j31);
fprintf("jacobian[3][2] = %s;\n", j32);
fprintf("jacobian[3][3] = %s;\n", j33);
fprintf("jacobian[3][4] = %s;\n", j34);
fprintf("jacobian[3][5] = %s;\n", j35);
fprintf("jacobian[4][0] = %s;\n", j40);
fprintf("jacobian[4][1] = %s;\n", j41);
fprintf("jacobian[4][2] = %s;\n", j42);
fprintf("jacobian[4][3] = %s;\n", j43);
fprintf("jacobian[4][4] = %s;\n", j44);
fprintf("jacobian[4][5] = %s;\n", j45);
fprintf("jacobian[5][0] = %s;\n", j50);
fprintf("jacobian[5][1] = %s;\n", j51);
fprintf("jacobian[5][2] = %s;\n", j52);
fprintf("jacobian[5][3] = %s;\n", j53);
fprintf("jacobian[5][4] = %s;\n", j54);
fprintf("jacobian[5][5] = %s;\n", j55);
