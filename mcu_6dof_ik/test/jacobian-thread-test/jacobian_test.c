#include <stdio.h>
//#include "linalg.h"
#include <math.h>
// gcc jacobian_test.c -lm

#define l1 1
#define l2 1
#define l3 1
#define l4 1
#define l5 1
#define l6 1
typedef long double type_t;

void print_jacobian_66(type_t mat[6][6])
{
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            printf("%5.2Lf ", mat[i][j]);
        }
        printf("\n");
    }
}

// int main(void)
// {
//     //get_jacobian_with_pthread(jacobian, cur_angle);
//     get_jacobian(jacobian, cur_angle);
//     return 0;
// }

int main(void)
{
    type_t jacobian[6][6] = {0};
    type_t cur_angle[6] = {1.4, 1, 0.3, 2, 0.3, 1};
    int is_out = 0;
    while (1)
    {
        printf("exit?(-1): ");
        scanf("%d", &is_out);
        if (is_out == -1)
            break;

        printf("angle (a0 a1...): ");
        scanf("%Lf %Lf %Lf %Lf %Lf %Lf", &cur_angle[0], &cur_angle[1], &cur_angle[2], &cur_angle[3], &cur_angle[4], &cur_angle[5]);

        type_t cos_a0 = cos(cur_angle[0]);
        type_t cos_a1 = cos(cur_angle[1]);
        type_t cos_a2 = cos(cur_angle[2]);
        type_t cos_a3 = cos(cur_angle[3]);
        type_t cos_a4 = cos(cur_angle[4]);
        type_t cos_a5 = cos(cur_angle[5]);

        type_t sin_a0 = sin(cur_angle[0]);
        type_t sin_a1 = sin(cur_angle[1]);
        type_t sin_a2 = sin(cur_angle[2]);
        type_t sin_a3 = sin(cur_angle[3]);
        type_t sin_a4 = sin(cur_angle[4]);
        type_t sin_a5 = sin(cur_angle[5]);

        jacobian[0][0] = -l5 * (cos_a4 * (cos_a0 * cos_a2 - cos_a1 * sin_a0 * sin_a2) - sin_a4 * (cos_a3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) - sin_a0 * sin_a1 * sin_a3)) - l6 * (sin_a5 * (sin_a4 * (cos_a0 * cos_a2 - cos_a1 * sin_a0 * sin_a2) + cos_a4 * (cos_a3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) - sin_a0 * sin_a1 * sin_a3)) + cos_a5 * (sin_a3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) + cos_a3 * sin_a0 * sin_a1)) - l4 * (sin_a3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) + cos_a3 * sin_a0 * sin_a1) - l1 * cos_a0 - l3 * (cos_a0 * cos_a2 - cos_a1 * sin_a0 * sin_a2) - l2 * sin_a0 * sin_a1;
        jacobian[0][1] = l4 * (cos_a0 * cos_a1 * cos_a3 - cos_a0 * cos_a2 * sin_a1 * sin_a3) + l5 * (sin_a4 * (cos_a0 * cos_a1 * sin_a3 + cos_a0 * cos_a2 * cos_a3 * sin_a1) + cos_a0 * cos_a4 * sin_a1 * sin_a2) + l6 * (cos_a5 * (cos_a0 * cos_a1 * cos_a3 - cos_a0 * cos_a2 * sin_a1 * sin_a3) - sin_a5 * (cos_a4 * (cos_a0 * cos_a1 * sin_a3 + cos_a0 * cos_a2 * cos_a3 * sin_a1) - cos_a0 * sin_a1 * sin_a2 * sin_a4)) + l2 * cos_a0 * cos_a1 + l3 * cos_a0 * sin_a1 * sin_a2;
        jacobian[0][2] = l6 * (sin_a5 * (sin_a4 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) - cos_a3 * cos_a4 * (cos_a2 * sin_a0 + cos_a0 * cos_a1 * sin_a2)) - cos_a5 * sin_a3 * (cos_a2 * sin_a0 + cos_a0 * cos_a1 * sin_a2)) + l5 * (cos_a4 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) + cos_a3 * sin_a4 * (cos_a2 * sin_a0 + cos_a0 * cos_a1 * sin_a2)) + l3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) - l4 * sin_a3 * (cos_a2 * sin_a0 + cos_a0 * cos_a1 * sin_a2);
        jacobian[0][3] = -l6 * (cos_a5 * (cos_a3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) + cos_a0 * sin_a1 * sin_a3) - cos_a4 * sin_a5 * (sin_a3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) - cos_a0 * cos_a3 * sin_a1)) - l4 * (cos_a3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) + cos_a0 * sin_a1 * sin_a3) - l5 * sin_a4 * (sin_a3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) - cos_a0 * cos_a3 * sin_a1);
        jacobian[0][4] = l5 * (sin_a4 * (cos_a2 * sin_a0 + cos_a0 * cos_a1 * sin_a2) + cos_a4 * (cos_a3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) + cos_a0 * sin_a1 * sin_a3)) - l6 * sin_a5 * (cos_a4 * (cos_a2 * sin_a0 + cos_a0 * cos_a1 * sin_a2) - sin_a4 * (cos_a3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) + cos_a0 * sin_a1 * sin_a3));
        jacobian[0][5] = -l6 * (cos_a5 * (sin_a4 * (cos_a2 * sin_a0 + cos_a0 * cos_a1 * sin_a2) + cos_a4 * (cos_a3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) + cos_a0 * sin_a1 * sin_a3)) - sin_a5 * (sin_a3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) - cos_a0 * cos_a3 * sin_a1));
        jacobian[1][0] = l2 * cos_a0 * sin_a1 - l4 * (sin_a3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) - cos_a0 * cos_a3 * sin_a1) - l1 * sin_a0 - l3 * (cos_a2 * sin_a0 + cos_a0 * cos_a1 * sin_a2) - l5 * (cos_a4 * (cos_a2 * sin_a0 + cos_a0 * cos_a1 * sin_a2) - sin_a4 * (cos_a3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) + cos_a0 * sin_a1 * sin_a3)) - l6 * (sin_a5 * (sin_a4 * (cos_a2 * sin_a0 + cos_a0 * cos_a1 * sin_a2) + cos_a4 * (cos_a3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) + cos_a0 * sin_a1 * sin_a3)) + cos_a5 * (sin_a3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) - cos_a0 * cos_a3 * sin_a1));
        jacobian[1][1] = l4 * (cos_a1 * cos_a3 * sin_a0 - cos_a2 * sin_a0 * sin_a1 * sin_a3) - l6 * (sin_a5 * (cos_a4 * (cos_a1 * sin_a0 * sin_a3 + cos_a2 * cos_a3 * sin_a0 * sin_a1) - sin_a0 * sin_a1 * sin_a2 * sin_a4) - cos_a5 * (cos_a1 * cos_a3 * sin_a0 - cos_a2 * sin_a0 * sin_a1 * sin_a3)) + l5 * (sin_a4 * (cos_a1 * sin_a0 * sin_a3 + cos_a2 * cos_a3 * sin_a0 * sin_a1) + cos_a4 * sin_a0 * sin_a1 * sin_a2) + l2 * cos_a1 * sin_a0 + l3 * sin_a0 * sin_a1 * sin_a2;
        jacobian[1][2] = l4 * sin_a3 * (cos_a0 * cos_a2 - cos_a1 * sin_a0 * sin_a2) - l5 * (cos_a4 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) + cos_a3 * sin_a4 * (cos_a0 * cos_a2 - cos_a1 * sin_a0 * sin_a2)) - l3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) - l6 * (sin_a5 * (sin_a4 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) - cos_a3 * cos_a4 * (cos_a0 * cos_a2 - cos_a1 * sin_a0 * sin_a2)) - cos_a5 * sin_a3 * (cos_a0 * cos_a2 - cos_a1 * sin_a0 * sin_a2));
        jacobian[1][3] = l6 * (cos_a5 * (cos_a3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) - sin_a0 * sin_a1 * sin_a3) - cos_a4 * sin_a5 * (sin_a3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) + cos_a3 * sin_a0 * sin_a1)) + l4 * (cos_a3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) - sin_a0 * sin_a1 * sin_a3) + l5 * sin_a4 * (sin_a3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) + cos_a3 * sin_a0 * sin_a1);
        jacobian[1][4] = l6 * sin_a5 * (cos_a4 * (cos_a0 * cos_a2 - cos_a1 * sin_a0 * sin_a2) - sin_a4 * (cos_a3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) - sin_a0 * sin_a1 * sin_a3)) - l5 * (sin_a4 * (cos_a0 * cos_a2 - cos_a1 * sin_a0 * sin_a2) + cos_a4 * (cos_a3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) - sin_a0 * sin_a1 * sin_a3));
        jacobian[1][5] = l6 * (cos_a5 * (sin_a4 * (cos_a0 * cos_a2 - cos_a1 * sin_a0 * sin_a2) + cos_a4 * (cos_a3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) - sin_a0 * sin_a1 * sin_a3)) - sin_a5 * (sin_a3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) + cos_a3 * sin_a0 * sin_a1));
        jacobian[2][0] = 0;
        jacobian[2][1] = l3 * cos_a1 * sin_a2 - l2 * sin_a1 - l4 * (cos_a3 * sin_a1 + cos_a1 * cos_a2 * sin_a3) - l6 * (cos_a5 * (cos_a3 * sin_a1 + cos_a1 * cos_a2 * sin_a3) - sin_a5 * (cos_a4 * (sin_a1 * sin_a3 - cos_a1 * cos_a2 * cos_a3) + cos_a1 * sin_a2 * sin_a4)) - l5 * (sin_a4 * (sin_a1 * sin_a3 - cos_a1 * cos_a2 * cos_a3) - cos_a1 * cos_a4 * sin_a2);
        jacobian[2][2] = l5 * (cos_a2 * cos_a4 * sin_a1 - cos_a3 * sin_a1 * sin_a2 * sin_a4) + l6 * (sin_a5 * (cos_a2 * sin_a1 * sin_a4 + cos_a3 * cos_a4 * sin_a1 * sin_a2) + cos_a5 * sin_a1 * sin_a2 * sin_a3) + l3 * cos_a2 * sin_a1 + l4 * sin_a1 * sin_a2 * sin_a3;
        jacobian[2][3] = l5 * sin_a4 * (cos_a1 * cos_a3 - cos_a2 * sin_a1 * sin_a3) - l4 * (cos_a1 * sin_a3 + cos_a2 * cos_a3 * sin_a1) - l6 * (cos_a5 * (cos_a1 * sin_a3 + cos_a2 * cos_a3 * sin_a1) + cos_a4 * sin_a5 * (cos_a1 * cos_a3 - cos_a2 * sin_a1 * sin_a3));
        jacobian[2][4] = l5 * (cos_a4 * (cos_a1 * sin_a3 + cos_a2 * cos_a3 * sin_a1) - sin_a1 * sin_a2 * sin_a4) + l6 * sin_a5 * (sin_a4 * (cos_a1 * sin_a3 + cos_a2 * cos_a3 * sin_a1) + cos_a4 * sin_a1 * sin_a2);
        jacobian[2][5] = -l6 * (sin_a5 * (cos_a1 * cos_a3 - cos_a2 * sin_a1 * sin_a3) + cos_a5 * (cos_a4 * (cos_a1 * sin_a3 + cos_a2 * cos_a3 * sin_a1) - sin_a1 * sin_a2 * sin_a4));
        jacobian[3][0] = -sin_a5 * (sin_a4 * (cos_a0 * cos_a2 - cos_a1 * sin_a0 * sin_a2) + cos_a4 * (cos_a3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) - sin_a0 * sin_a1 * sin_a3)) - cos_a5 * (sin_a3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) + cos_a3 * sin_a0 * sin_a1);
        jacobian[3][1] = cos_a5 * (cos_a0 * cos_a1 * cos_a3 - cos_a0 * cos_a2 * sin_a1 * sin_a3) - sin_a5 * (cos_a4 * (cos_a0 * cos_a1 * sin_a3 + cos_a0 * cos_a2 * cos_a3 * sin_a1) - cos_a0 * sin_a1 * sin_a2 * sin_a4);
        jacobian[3][2] = sin_a5 * (sin_a4 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) - cos_a3 * cos_a4 * (cos_a2 * sin_a0 + cos_a0 * cos_a1 * sin_a2)) - cos_a5 * sin_a3 * (cos_a2 * sin_a0 + cos_a0 * cos_a1 * sin_a2);
        jacobian[3][3] = cos_a4 * sin_a5 * (sin_a3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) - cos_a0 * cos_a3 * sin_a1) - cos_a5 * (cos_a3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) + cos_a0 * sin_a1 * sin_a3);
        jacobian[3][4] = -sin_a5 * (cos_a4 * (cos_a2 * sin_a0 + cos_a0 * cos_a1 * sin_a2) - sin_a4 * (cos_a3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) + cos_a0 * sin_a1 * sin_a3));
        jacobian[3][5] = sin_a5 * (sin_a3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) - cos_a0 * cos_a3 * sin_a1) - cos_a5 * (sin_a4 * (cos_a2 * sin_a0 + cos_a0 * cos_a1 * sin_a2) + cos_a4 * (cos_a3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) + cos_a0 * sin_a1 * sin_a3));
        jacobian[4][0] = -sin_a5 * (sin_a4 * (cos_a2 * sin_a0 + cos_a0 * cos_a1 * sin_a2) + cos_a4 * (cos_a3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) + cos_a0 * sin_a1 * sin_a3)) - cos_a5 * (sin_a3 * (sin_a0 * sin_a2 - cos_a0 * cos_a1 * cos_a2) - cos_a0 * cos_a3 * sin_a1);
        jacobian[4][1] = cos_a5 * (cos_a1 * cos_a3 * sin_a0 - cos_a2 * sin_a0 * sin_a1 * sin_a3) - sin_a5 * (cos_a4 * (cos_a1 * sin_a0 * sin_a3 + cos_a2 * cos_a3 * sin_a0 * sin_a1) - sin_a0 * sin_a1 * sin_a2 * sin_a4);
        jacobian[4][2] = cos_a5 * sin_a3 * (cos_a0 * cos_a2 - cos_a1 * sin_a0 * sin_a2) - sin_a5 * (sin_a4 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) - cos_a3 * cos_a4 * (cos_a0 * cos_a2 - cos_a1 * sin_a0 * sin_a2));
        jacobian[4][3] = cos_a5 * (cos_a3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) - sin_a0 * sin_a1 * sin_a3) - cos_a4 * sin_a5 * (sin_a3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) + cos_a3 * sin_a0 * sin_a1);
        jacobian[4][4] = sin_a5 * (cos_a4 * (cos_a0 * cos_a2 - cos_a1 * sin_a0 * sin_a2) - sin_a4 * (cos_a3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) - sin_a0 * sin_a1 * sin_a3));
        jacobian[4][5] = cos_a5 * (sin_a4 * (cos_a0 * cos_a2 - cos_a1 * sin_a0 * sin_a2) + cos_a4 * (cos_a3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) - sin_a0 * sin_a1 * sin_a3)) - sin_a5 * (sin_a3 * (cos_a0 * sin_a2 + cos_a1 * cos_a2 * sin_a0) + cos_a3 * sin_a0 * sin_a1);
        jacobian[5][0] = 0;
        jacobian[5][1] = sin_a5 * (cos_a4 * (sin_a1 * sin_a3 - cos_a1 * cos_a2 * cos_a3) + cos_a1 * sin_a2 * sin_a4) - cos_a5 * (cos_a3 * sin_a1 + cos_a1 * cos_a2 * sin_a3);
        jacobian[5][2] = sin_a5 * (cos_a2 * sin_a1 * sin_a4 + cos_a3 * cos_a4 * sin_a1 * sin_a2) + cos_a5 * sin_a1 * sin_a2 * sin_a3;
        jacobian[5][3] = -cos_a5 * (cos_a1 * sin_a3 + cos_a2 * cos_a3 * sin_a1) - cos_a4 * sin_a5 * (cos_a1 * cos_a3 - cos_a2 * sin_a1 * sin_a3);
        jacobian[5][4] = sin_a5 * (sin_a4 * (cos_a1 * sin_a3 + cos_a2 * cos_a3 * sin_a1) + cos_a4 * sin_a1 * sin_a2);
        jacobian[5][5] = -sin_a5 * (cos_a1 * cos_a3 - cos_a2 * sin_a1 * sin_a3) - cos_a5 * (cos_a4 * (cos_a1 * sin_a3 + cos_a2 * cos_a3 * sin_a1) - sin_a1 * sin_a2 * sin_a4);
        print_jacobian_66(jacobian);
    }
    return 0;
}
