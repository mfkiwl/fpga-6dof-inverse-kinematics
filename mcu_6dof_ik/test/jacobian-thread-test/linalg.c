#include "linalg.h"

int g_dim = 0;

type_t (*g_mat_a)[MAX_LEN] = NULL;
type_t (*g_mat_b)[MAX_LEN] = NULL;
type_t (*g_mat_c)[MAX_LEN] = NULL;
type_t g_trans_mat_b[MAX_LEN][MAX_LEN] = {0};

type_t g_cos_a0 = 0, g_cos_a1 = 0, g_cos_a2 = 0, g_cos_a3 = 0, g_cos_a4 = 0, g_cos_a5 = 0;
type_t g_sin_a0 = 0, g_sin_a1 = 0, g_sin_a2 = 0, g_sin_a3 = 0, g_sin_a4 = 0, g_sin_a5 = 0;
type_t **g_jacobian = NULL;
bool did_init_g_get_jacobian_func_arr = true;
void* (*g_get_jacobian_func_arr[6][6])(void*) = {NULL};

void *get_jacobian_00(void* arg)
{
    ((type_t(*)[6]) arg)[0][0] = -l5 * (g_cos_a4 * (g_cos_a0 * g_cos_a2 - g_cos_a1 * g_sin_a0 * g_sin_a2) - g_sin_a4 * (g_cos_a3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) - g_sin_a0 * g_sin_a1 * g_sin_a3)) - l6 * (g_sin_a5 * (g_sin_a4 * (g_cos_a0 * g_cos_a2 - g_cos_a1 * g_sin_a0 * g_sin_a2) + g_cos_a4 * (g_cos_a3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) - g_sin_a0 * g_sin_a1 * g_sin_a3)) + g_cos_a5 * (g_sin_a3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) + g_cos_a3 * g_sin_a0 * g_sin_a1)) - l4 * (g_sin_a3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) + g_cos_a3 * g_sin_a0 * g_sin_a1) - l1 * g_cos_a0 - l3 * (g_cos_a0 * g_cos_a2 - g_cos_a1 * g_sin_a0 * g_sin_a2) - l2 * g_sin_a0 * g_sin_a1;
    pthread_exit(NULL);
}
void *get_jacobian_01(void* arg)
{
    ((type_t(*)[6]) arg)[0][1] = l4 * (g_cos_a0 * g_cos_a1 * g_cos_a3 - g_cos_a0 * g_cos_a2 * g_sin_a1 * g_sin_a3) + l5 * (g_sin_a4 * (g_cos_a0 * g_cos_a1 * g_sin_a3 + g_cos_a0 * g_cos_a2 * g_cos_a3 * g_sin_a1) + g_cos_a0 * g_cos_a4 * g_sin_a1 * g_sin_a2) + l6 * (g_cos_a5 * (g_cos_a0 * g_cos_a1 * g_cos_a3 - g_cos_a0 * g_cos_a2 * g_sin_a1 * g_sin_a3) - g_sin_a5 * (g_cos_a4 * (g_cos_a0 * g_cos_a1 * g_sin_a3 + g_cos_a0 * g_cos_a2 * g_cos_a3 * g_sin_a1) - g_cos_a0 * g_sin_a1 * g_sin_a2 * g_sin_a4)) + l2 * g_cos_a0 * g_cos_a1 + l3 * g_cos_a0 * g_sin_a1 * g_sin_a2;
    pthread_exit(NULL);
}
void *get_jacobian_02(void* arg)
{
    ((type_t(*)[6]) arg)[0][2] = l6 * (g_sin_a5 * (g_sin_a4 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) - g_cos_a3 * g_cos_a4 * (g_cos_a2 * g_sin_a0 + g_cos_a0 * g_cos_a1 * g_sin_a2)) - g_cos_a5 * g_sin_a3 * (g_cos_a2 * g_sin_a0 + g_cos_a0 * g_cos_a1 * g_sin_a2)) + l5 * (g_cos_a4 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) + g_cos_a3 * g_sin_a4 * (g_cos_a2 * g_sin_a0 + g_cos_a0 * g_cos_a1 * g_sin_a2)) + l3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) - l4 * g_sin_a3 * (g_cos_a2 * g_sin_a0 + g_cos_a0 * g_cos_a1 * g_sin_a2);
    pthread_exit(NULL);
}
void *get_jacobian_03(void* arg)
{
    ((type_t(*)[6]) arg)[0][3] = -l6 * (g_cos_a5 * (g_cos_a3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) + g_cos_a0 * g_sin_a1 * g_sin_a3) - g_cos_a4 * g_sin_a5 * (g_sin_a3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) - g_cos_a0 * g_cos_a3 * g_sin_a1)) - l4 * (g_cos_a3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) + g_cos_a0 * g_sin_a1 * g_sin_a3) - l5 * g_sin_a4 * (g_sin_a3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) - g_cos_a0 * g_cos_a3 * g_sin_a1);
    pthread_exit(NULL);
}
void *get_jacobian_04(void* arg)
{
    ((type_t(*)[6]) arg)[0][4] = l5 * (g_sin_a4 * (g_cos_a2 * g_sin_a0 + g_cos_a0 * g_cos_a1 * g_sin_a2) + g_cos_a4 * (g_cos_a3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) + g_cos_a0 * g_sin_a1 * g_sin_a3)) - l6 * g_sin_a5 * (g_cos_a4 * (g_cos_a2 * g_sin_a0 + g_cos_a0 * g_cos_a1 * g_sin_a2) - g_sin_a4 * (g_cos_a3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) + g_cos_a0 * g_sin_a1 * g_sin_a3));
    pthread_exit(NULL);
}
void *get_jacobian_05(void* arg)
{
    ((type_t(*)[6]) arg)[0][5] = -l6 * (g_cos_a5 * (g_sin_a4 * (g_cos_a2 * g_sin_a0 + g_cos_a0 * g_cos_a1 * g_sin_a2) + g_cos_a4 * (g_cos_a3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) + g_cos_a0 * g_sin_a1 * g_sin_a3)) - g_sin_a5 * (g_sin_a3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) - g_cos_a0 * g_cos_a3 * g_sin_a1));
    pthread_exit(NULL);
}
void *get_jacobian_10(void* arg)
{
    ((type_t(*)[6]) arg)[1][0] = l2 * g_cos_a0 * g_sin_a1 - l4 * (g_sin_a3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) - g_cos_a0 * g_cos_a3 * g_sin_a1) - l1 * g_sin_a0 - l3 * (g_cos_a2 * g_sin_a0 + g_cos_a0 * g_cos_a1 * g_sin_a2) - l5 * (g_cos_a4 * (g_cos_a2 * g_sin_a0 + g_cos_a0 * g_cos_a1 * g_sin_a2) - g_sin_a4 * (g_cos_a3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) + g_cos_a0 * g_sin_a1 * g_sin_a3)) - l6 * (g_sin_a5 * (g_sin_a4 * (g_cos_a2 * g_sin_a0 + g_cos_a0 * g_cos_a1 * g_sin_a2) + g_cos_a4 * (g_cos_a3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) + g_cos_a0 * g_sin_a1 * g_sin_a3)) + g_cos_a5 * (g_sin_a3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) - g_cos_a0 * g_cos_a3 * g_sin_a1));
    pthread_exit(NULL);
}
void *get_jacobian_11(void* arg)
{
    ((type_t(*)[6]) arg)[1][1] = l4 * (g_cos_a1 * g_cos_a3 * g_sin_a0 - g_cos_a2 * g_sin_a0 * g_sin_a1 * g_sin_a3) - l6 * (g_sin_a5 * (g_cos_a4 * (g_cos_a1 * g_sin_a0 * g_sin_a3 + g_cos_a2 * g_cos_a3 * g_sin_a0 * g_sin_a1) - g_sin_a0 * g_sin_a1 * g_sin_a2 * g_sin_a4) - g_cos_a5 * (g_cos_a1 * g_cos_a3 * g_sin_a0 - g_cos_a2 * g_sin_a0 * g_sin_a1 * g_sin_a3)) + l5 * (g_sin_a4 * (g_cos_a1 * g_sin_a0 * g_sin_a3 + g_cos_a2 * g_cos_a3 * g_sin_a0 * g_sin_a1) + g_cos_a4 * g_sin_a0 * g_sin_a1 * g_sin_a2) + l2 * g_cos_a1 * g_sin_a0 + l3 * g_sin_a0 * g_sin_a1 * g_sin_a2;
    pthread_exit(NULL);
}
void *get_jacobian_12(void* arg)
{
    ((type_t(*)[6]) arg)[1][2] = l4 * g_sin_a3 * (g_cos_a0 * g_cos_a2 - g_cos_a1 * g_sin_a0 * g_sin_a2) - l5 * (g_cos_a4 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) + g_cos_a3 * g_sin_a4 * (g_cos_a0 * g_cos_a2 - g_cos_a1 * g_sin_a0 * g_sin_a2)) - l3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) - l6 * (g_sin_a5 * (g_sin_a4 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) - g_cos_a3 * g_cos_a4 * (g_cos_a0 * g_cos_a2 - g_cos_a1 * g_sin_a0 * g_sin_a2)) - g_cos_a5 * g_sin_a3 * (g_cos_a0 * g_cos_a2 - g_cos_a1 * g_sin_a0 * g_sin_a2));
    pthread_exit(NULL);
}
void *get_jacobian_13(void* arg)
{
    ((type_t(*)[6]) arg)[1][3] = l6 * (g_cos_a5 * (g_cos_a3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) - g_sin_a0 * g_sin_a1 * g_sin_a3) - g_cos_a4 * g_sin_a5 * (g_sin_a3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) + g_cos_a3 * g_sin_a0 * g_sin_a1)) + l4 * (g_cos_a3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) - g_sin_a0 * g_sin_a1 * g_sin_a3) + l5 * g_sin_a4 * (g_sin_a3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) + g_cos_a3 * g_sin_a0 * g_sin_a1);
    pthread_exit(NULL);
}
void *get_jacobian_14(void* arg)
{
    ((type_t(*)[6]) arg)[1][4] = l6 * g_sin_a5 * (g_cos_a4 * (g_cos_a0 * g_cos_a2 - g_cos_a1 * g_sin_a0 * g_sin_a2) - g_sin_a4 * (g_cos_a3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) - g_sin_a0 * g_sin_a1 * g_sin_a3)) - l5 * (g_sin_a4 * (g_cos_a0 * g_cos_a2 - g_cos_a1 * g_sin_a0 * g_sin_a2) + g_cos_a4 * (g_cos_a3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) - g_sin_a0 * g_sin_a1 * g_sin_a3));
    pthread_exit(NULL);
}
void *get_jacobian_15(void* arg)
{
    ((type_t(*)[6]) arg)[1][5] = l6 * (g_cos_a5 * (g_sin_a4 * (g_cos_a0 * g_cos_a2 - g_cos_a1 * g_sin_a0 * g_sin_a2) + g_cos_a4 * (g_cos_a3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) - g_sin_a0 * g_sin_a1 * g_sin_a3)) - g_sin_a5 * (g_sin_a3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) + g_cos_a3 * g_sin_a0 * g_sin_a1));
    pthread_exit(NULL);
}
void *get_jacobian_20(void* arg)
{
    ((type_t(*)[6]) arg)[2][0] = 0;
    pthread_exit(NULL);
}
void *get_jacobian_21(void* arg)
{
    ((type_t(*)[6]) arg)[2][1] = l3 * g_cos_a1 * g_sin_a2 - l2 * g_sin_a1 - l4 * (g_cos_a3 * g_sin_a1 + g_cos_a1 * g_cos_a2 * g_sin_a3) - l6 * (g_cos_a5 * (g_cos_a3 * g_sin_a1 + g_cos_a1 * g_cos_a2 * g_sin_a3) - g_sin_a5 * (g_cos_a4 * (g_sin_a1 * g_sin_a3 - g_cos_a1 * g_cos_a2 * g_cos_a3) + g_cos_a1 * g_sin_a2 * g_sin_a4)) - l5 * (g_sin_a4 * (g_sin_a1 * g_sin_a3 - g_cos_a1 * g_cos_a2 * g_cos_a3) - g_cos_a1 * g_cos_a4 * g_sin_a2);
    pthread_exit(NULL);
}
void *get_jacobian_22(void* arg)
{
    ((type_t(*)[6]) arg)[2][2] = l5 * (g_cos_a2 * g_cos_a4 * g_sin_a1 - g_cos_a3 * g_sin_a1 * g_sin_a2 * g_sin_a4) + l6 * (g_sin_a5 * (g_cos_a2 * g_sin_a1 * g_sin_a4 + g_cos_a3 * g_cos_a4 * g_sin_a1 * g_sin_a2) + g_cos_a5 * g_sin_a1 * g_sin_a2 * g_sin_a3) + l3 * g_cos_a2 * g_sin_a1 + l4 * g_sin_a1 * g_sin_a2 * g_sin_a3;
    pthread_exit(NULL);
}
void *get_jacobian_23(void* arg)
{
    ((type_t(*)[6]) arg)[2][3] = l5 * g_sin_a4 * (g_cos_a1 * g_cos_a3 - g_cos_a2 * g_sin_a1 * g_sin_a3) - l4 * (g_cos_a1 * g_sin_a3 + g_cos_a2 * g_cos_a3 * g_sin_a1) - l6 * (g_cos_a5 * (g_cos_a1 * g_sin_a3 + g_cos_a2 * g_cos_a3 * g_sin_a1) + g_cos_a4 * g_sin_a5 * (g_cos_a1 * g_cos_a3 - g_cos_a2 * g_sin_a1 * g_sin_a3));
    pthread_exit(NULL);
}
void *get_jacobian_24(void* arg)
{
    ((type_t(*)[6]) arg)[2][4] = l5 * (g_cos_a4 * (g_cos_a1 * g_sin_a3 + g_cos_a2 * g_cos_a3 * g_sin_a1) - g_sin_a1 * g_sin_a2 * g_sin_a4) + l6 * g_sin_a5 * (g_sin_a4 * (g_cos_a1 * g_sin_a3 + g_cos_a2 * g_cos_a3 * g_sin_a1) + g_cos_a4 * g_sin_a1 * g_sin_a2);
    pthread_exit(NULL);
}
void *get_jacobian_25(void* arg)
{
    ((type_t(*)[6]) arg)[2][5] = -l6 * (g_sin_a5 * (g_cos_a1 * g_cos_a3 - g_cos_a2 * g_sin_a1 * g_sin_a3) + g_cos_a5 * (g_cos_a4 * (g_cos_a1 * g_sin_a3 + g_cos_a2 * g_cos_a3 * g_sin_a1) - g_sin_a1 * g_sin_a2 * g_sin_a4));
    pthread_exit(NULL);
}
void *get_jacobian_30(void* arg)
{
    ((type_t(*)[6]) arg)[3][0] = -g_sin_a5 * (g_sin_a4 * (g_cos_a0 * g_cos_a2 - g_cos_a1 * g_sin_a0 * g_sin_a2) + g_cos_a4 * (g_cos_a3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) - g_sin_a0 * g_sin_a1 * g_sin_a3)) - g_cos_a5 * (g_sin_a3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) + g_cos_a3 * g_sin_a0 * g_sin_a1);
    pthread_exit(NULL);
}
void *get_jacobian_31(void* arg)
{
    ((type_t(*)[6]) arg)[3][1] = g_cos_a5 * (g_cos_a0 * g_cos_a1 * g_cos_a3 - g_cos_a0 * g_cos_a2 * g_sin_a1 * g_sin_a3) - g_sin_a5 * (g_cos_a4 * (g_cos_a0 * g_cos_a1 * g_sin_a3 + g_cos_a0 * g_cos_a2 * g_cos_a3 * g_sin_a1) - g_cos_a0 * g_sin_a1 * g_sin_a2 * g_sin_a4);
    pthread_exit(NULL);
}
void *get_jacobian_32(void* arg)
{
    ((type_t(*)[6]) arg)[3][2] = g_sin_a5 * (g_sin_a4 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) - g_cos_a3 * g_cos_a4 * (g_cos_a2 * g_sin_a0 + g_cos_a0 * g_cos_a1 * g_sin_a2)) - g_cos_a5 * g_sin_a3 * (g_cos_a2 * g_sin_a0 + g_cos_a0 * g_cos_a1 * g_sin_a2);
    pthread_exit(NULL);
}
void *get_jacobian_33(void* arg)
{
    ((type_t(*)[6]) arg)[3][3] = g_cos_a4 * g_sin_a5 * (g_sin_a3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) - g_cos_a0 * g_cos_a3 * g_sin_a1) - g_cos_a5 * (g_cos_a3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) + g_cos_a0 * g_sin_a1 * g_sin_a3);
    pthread_exit(NULL);
}
void *get_jacobian_34(void* arg)
{
    ((type_t(*)[6]) arg)[3][4] = -g_sin_a5 * (g_cos_a4 * (g_cos_a2 * g_sin_a0 + g_cos_a0 * g_cos_a1 * g_sin_a2) - g_sin_a4 * (g_cos_a3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) + g_cos_a0 * g_sin_a1 * g_sin_a3));
    pthread_exit(NULL);
}
void *get_jacobian_35(void* arg)
{
    ((type_t(*)[6]) arg)[3][5] = g_sin_a5 * (g_sin_a3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) - g_cos_a0 * g_cos_a3 * g_sin_a1) - g_cos_a5 * (g_sin_a4 * (g_cos_a2 * g_sin_a0 + g_cos_a0 * g_cos_a1 * g_sin_a2) + g_cos_a4 * (g_cos_a3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) + g_cos_a0 * g_sin_a1 * g_sin_a3));
    pthread_exit(NULL);
}
void *get_jacobian_40(void* arg)
{
    ((type_t(*)[6]) arg)[4][0] = -g_sin_a5 * (g_sin_a4 * (g_cos_a2 * g_sin_a0 + g_cos_a0 * g_cos_a1 * g_sin_a2) + g_cos_a4 * (g_cos_a3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) + g_cos_a0 * g_sin_a1 * g_sin_a3)) - g_cos_a5 * (g_sin_a3 * (g_sin_a0 * g_sin_a2 - g_cos_a0 * g_cos_a1 * g_cos_a2) - g_cos_a0 * g_cos_a3 * g_sin_a1);
    pthread_exit(NULL);
}
void *get_jacobian_41(void* arg)
{
    ((type_t(*)[6]) arg)[4][1] = g_cos_a5 * (g_cos_a1 * g_cos_a3 * g_sin_a0 - g_cos_a2 * g_sin_a0 * g_sin_a1 * g_sin_a3) - g_sin_a5 * (g_cos_a4 * (g_cos_a1 * g_sin_a0 * g_sin_a3 + g_cos_a2 * g_cos_a3 * g_sin_a0 * g_sin_a1) - g_sin_a0 * g_sin_a1 * g_sin_a2 * g_sin_a4);
    pthread_exit(NULL);
}
void *get_jacobian_42(void* arg)
{
    ((type_t(*)[6]) arg)[4][2] = g_cos_a5 * g_sin_a3 * (g_cos_a0 * g_cos_a2 - g_cos_a1 * g_sin_a0 * g_sin_a2) - g_sin_a5 * (g_sin_a4 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) - g_cos_a3 * g_cos_a4 * (g_cos_a0 * g_cos_a2 - g_cos_a1 * g_sin_a0 * g_sin_a2));
    pthread_exit(NULL);
}
void *get_jacobian_43(void* arg)
{
    ((type_t(*)[6]) arg)[4][3] = g_cos_a5 * (g_cos_a3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) - g_sin_a0 * g_sin_a1 * g_sin_a3) - g_cos_a4 * g_sin_a5 * (g_sin_a3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) + g_cos_a3 * g_sin_a0 * g_sin_a1);
    pthread_exit(NULL);
}
void *get_jacobian_44(void* arg)
{
    ((type_t(*)[6]) arg)[4][4] = g_sin_a5 * (g_cos_a4 * (g_cos_a0 * g_cos_a2 - g_cos_a1 * g_sin_a0 * g_sin_a2) - g_sin_a4 * (g_cos_a3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) - g_sin_a0 * g_sin_a1 * g_sin_a3));
    pthread_exit(NULL);
}
void *get_jacobian_45(void* arg)
{
    ((type_t(*)[6]) arg)[4][5] = g_cos_a5 * (g_sin_a4 * (g_cos_a0 * g_cos_a2 - g_cos_a1 * g_sin_a0 * g_sin_a2) + g_cos_a4 * (g_cos_a3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) - g_sin_a0 * g_sin_a1 * g_sin_a3)) - g_sin_a5 * (g_sin_a3 * (g_cos_a0 * g_sin_a2 + g_cos_a1 * g_cos_a2 * g_sin_a0) + g_cos_a3 * g_sin_a0 * g_sin_a1);
    pthread_exit(NULL);
}
void *get_jacobian_50(void* arg)
{
    ((type_t(*)[6]) arg)[5][0] = 0;
    pthread_exit(NULL);
}
void *get_jacobian_51(void* arg)
{
    ((type_t(*)[6]) arg)[5][1] = g_sin_a5 * (g_cos_a4 * (g_sin_a1 * g_sin_a3 - g_cos_a1 * g_cos_a2 * g_cos_a3) + g_cos_a1 * g_sin_a2 * g_sin_a4) - g_cos_a5 * (g_cos_a3 * g_sin_a1 + g_cos_a1 * g_cos_a2 * g_sin_a3);
    pthread_exit(NULL);
}
void *get_jacobian_52(void* arg)
{
    ((type_t(*)[6]) arg)[5][2] = g_sin_a5 * (g_cos_a2 * g_sin_a1 * g_sin_a4 + g_cos_a3 * g_cos_a4 * g_sin_a1 * g_sin_a2) + g_cos_a5 * g_sin_a1 * g_sin_a2 * g_sin_a3;
    pthread_exit(NULL);
}
void *get_jacobian_53(void* arg)
{
    ((type_t(*)[6]) arg)[5][3] = -g_cos_a5 * (g_cos_a1 * g_sin_a3 + g_cos_a2 * g_cos_a3 * g_sin_a1) - g_cos_a4 * g_sin_a5 * (g_cos_a1 * g_cos_a3 - g_cos_a2 * g_sin_a1 * g_sin_a3);
    pthread_exit(NULL);
}
void *get_jacobian_54(void* arg)
{
    ((type_t(*)[6]) arg)[5][4] = g_sin_a5 * (g_sin_a4 * (g_cos_a1 * g_sin_a3 + g_cos_a2 * g_cos_a3 * g_sin_a1) + g_cos_a4 * g_sin_a1 * g_sin_a2);
    pthread_exit(NULL);
}
void *get_jacobian_55(void* arg)
{
    ((type_t(*)[6]) arg)[5][5] = -g_sin_a5 * (g_cos_a1 * g_cos_a3 - g_cos_a2 * g_sin_a1 * g_sin_a3) - g_cos_a5 * (g_cos_a4 * (g_cos_a1 * g_sin_a3 + g_cos_a2 * g_cos_a3 * g_sin_a1) - g_sin_a1 * g_sin_a2 * g_sin_a4);
    pthread_exit(NULL);
}
void init_jacobian_func_arr(void)
{
    g_get_jacobian_func_arr[0][0] = get_jacobian_00;
    g_get_jacobian_func_arr[0][1] = get_jacobian_01;
    g_get_jacobian_func_arr[0][2] = get_jacobian_02;
    g_get_jacobian_func_arr[0][3] = get_jacobian_03;
    g_get_jacobian_func_arr[0][4] = get_jacobian_04;
    g_get_jacobian_func_arr[0][5] = get_jacobian_05;
    g_get_jacobian_func_arr[1][0] = get_jacobian_10;
    g_get_jacobian_func_arr[1][1] = get_jacobian_11;
    g_get_jacobian_func_arr[1][2] = get_jacobian_12;
    g_get_jacobian_func_arr[1][3] = get_jacobian_13;
    g_get_jacobian_func_arr[1][4] = get_jacobian_14;
    g_get_jacobian_func_arr[1][5] = get_jacobian_15;
    g_get_jacobian_func_arr[2][0] = get_jacobian_20;
    g_get_jacobian_func_arr[2][1] = get_jacobian_21;
    g_get_jacobian_func_arr[2][2] = get_jacobian_22;
    g_get_jacobian_func_arr[2][3] = get_jacobian_23;
    g_get_jacobian_func_arr[2][4] = get_jacobian_24;
    g_get_jacobian_func_arr[2][5] = get_jacobian_25;
    g_get_jacobian_func_arr[3][0] = get_jacobian_30;
    g_get_jacobian_func_arr[3][1] = get_jacobian_31;
    g_get_jacobian_func_arr[3][2] = get_jacobian_32;
    g_get_jacobian_func_arr[3][3] = get_jacobian_33;
    g_get_jacobian_func_arr[3][4] = get_jacobian_34;
    g_get_jacobian_func_arr[3][5] = get_jacobian_35;
    g_get_jacobian_func_arr[4][0] = get_jacobian_40;
    g_get_jacobian_func_arr[4][1] = get_jacobian_41;
    g_get_jacobian_func_arr[4][2] = get_jacobian_42;
    g_get_jacobian_func_arr[4][3] = get_jacobian_43;
    g_get_jacobian_func_arr[4][4] = get_jacobian_44;
    g_get_jacobian_func_arr[4][5] = get_jacobian_45;
    g_get_jacobian_func_arr[5][0] = get_jacobian_50;
    g_get_jacobian_func_arr[5][1] = get_jacobian_51;
    g_get_jacobian_func_arr[5][2] = get_jacobian_52;
    g_get_jacobian_func_arr[5][3] = get_jacobian_53;
    g_get_jacobian_func_arr[5][4] = get_jacobian_54;
    g_get_jacobian_func_arr[5][5] = get_jacobian_55;
}

void init_tigonometric_func(type_t cur_angle[6])
{
    g_cos_a0 = cos(cur_angle[0]);
    g_cos_a1 = cos(cur_angle[1]);
    g_cos_a2 = cos(cur_angle[2]);
    g_cos_a3 = cos(cur_angle[3]);
    g_cos_a4 = cos(cur_angle[4]);
    g_cos_a5 = cos(cur_angle[5]);

    g_sin_a0 = sin(cur_angle[0]);
    g_sin_a1 = sin(cur_angle[1]);
    g_sin_a2 = sin(cur_angle[2]);
    g_sin_a3 = sin(cur_angle[3]);
    g_sin_a4 = sin(cur_angle[4]);
    g_sin_a5 = sin(cur_angle[5]);
}

void get_jacobian_with_pthread(type_t jacobian[][6], type_t cur_angle[6])
{
    g_jacobian = (type_t**) jacobian;
    if (did_init_g_get_jacobian_func_arr)
    {
        init_jacobian_func_arr();
        did_init_g_get_jacobian_func_arr = false;
    }
    init_tigonometric_func(cur_angle);
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    pthread_t thread_pool[6][6];

    // assign B matrix value intotranposed B temporary matrix
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            pthread_create(&thread_pool[i][j], &attr, g_get_jacobian_func_arr[i][j], (void*) jacobian);
        }
    }

    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            pthread_join(thread_pool[i][j], NULL);
        }
    }
}

void get_jacobian(type_t jacobian[][6], type_t cur_angle[6])
{
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
}