// gcc -c linalg.c -lm -lpthread
// gcc jacobian_test.c linalg.o -lm -lpthread

#include <stdio.h>
#include "linalg.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#define PI 3.141592

void print_jacobian_66(type_t mat[6][6])
{
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            printf("%10lf\t", mat[i][j]);
        }
        printf("\n");
    }
}

int main(void)
{
    type_t jacobian[6][6] = {0};
    // angle in radian
    type_t cur_angle[6] = {0};
    char c;

    for (;;)
    {
        printf("exit ? (y/外) ");
        c = getchar();
        getchar();

        srand((unsigned)time(NULL));
        if (c == 'y')
            break;
        for (int i = 0; i < 6; i++)
        {
            cur_angle[i] = ((type_t) rand() / RAND_MAX) * PI;
            printf("a[%d]: %lf * pi radian\n", i, cur_angle[i]);
        }

        printf("without pthread\n");
        get_jacobian_with_pthread(jacobian, cur_angle);
        print_jacobian_66(jacobian);

        printf("\nwith pthread\n");
        get_jacobian(jacobian, cur_angle);
        print_jacobian_66(jacobian);
    }

    return 0;
}
