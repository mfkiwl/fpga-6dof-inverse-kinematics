#include "linalg.h"

int g_dim = 0;

type_t (*g_mat_a)[MAX_LEN] = NULL;
type_t (*g_mat_b)[MAX_LEN] = NULL;
type_t (*g_mat_c)[MAX_LEN] = NULL;
type_t g_trans_mat_b[MAX_LEN][MAX_LEN] = {0};

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

// ERROR CODE
//  1. DIV_BY_ZERO
//  2. DONE
int inverse_matrix(int dim, type_t inv_a[][MAX_LEN], type_t a[][MAX_LEN])
{
    type_t x[MAX_LEN], ratio;
    int i, j, k, n = dim;
    /*  Identity Matrix of Order n */
    memcpy(inv_a, a, sizeof(type_t) * MAX_LEN * MAX_LEN);

    for (i = 1; i <= n; i++)
    {
        for (j = 1; j <= n; j++)
        {
            if (i == j)
            {
                inv_a[i][j + n] = 1;
            }
            else
            {
                inv_a[i][j + n] = 0;
            }
        }
    }
    // Gauss Jordan Elimination
    for (i = 1; i <= n; i++)
    {
        if (a[i][i] == 0.0)
        {
            // error: div by zero
            return DIV_BY_ZERO;
        }
        for (j = 1; j <= n; j++)
        {
            if (i != j)
            {
                ratio = a[j][i] / a[i][i];
                for (k = 1; k <= 2 * n; k++)
                {
                    a[j][k] = a[j][k] - ratio * a[i][k];
                }
            }
        }
    }
    // Main diagonal to 1
    for (i = 1; i <= n; i++)
    {
        for (j = n + 1; j <= 2 * n; j++)
        {
            a[i][j] = a[i][j] / a[i][i];
        }
    }
    return DONE;
}

void *compute_row(void *args)
{
    arguments *argument = (arguments *)args;
    int start_row = argument->start_row;
    int end_row = argument->end_row;
    type_t sum = 0.0;
    for (int i = start_row; i < end_row + 1; i++)
    {
        for (int j = 0; j < g_dim; j++)
        {
            // OF
            sum = 0.0;
            for (int k = 0; k < g_dim; k++)
                sum += g_mat_a[i][k] * g_trans_mat_b[j][k];
            g_mat_c[i][j] = sum;
        }
    }
    pthread_exit(NULL);
}

void multiply(int dim, type_t mat_a[][MAX_LEN], type_t mat_b[][MAX_LEN], type_t mat_c[][MAX_LEN])
{
    g_dim = dim;
    g_mat_a = mat_a;
    g_mat_b = mat_b;
    g_mat_c = mat_c;

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    pthread_t thread_pool[NUM_T];
    arguments args_pool[NUM_T];

    // assign B matrix value intotranposed B temporary matrix
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < dim; j++)
        {
            g_trans_mat_b[j][i] = g_mat_b[i][j];
        }
    }

    int n_split = dim < NUM_T ? dim : NUM_T;
    int n_work = dim < NUM_T ? 1 : dim / NUM_T;

    for (int i = 0; i < n_split; i++)
    {
        arguments args;
        args.start_row = i * n_work;
        args.end_row = args.start_row + n_work;
        args_pool[i] = args;
        pthread_create(&thread_pool[i], &attr, compute_row, (void *)&args_pool[i]);
    }

    for (int i = 0; i < n_split; i++)
    {
        pthread_join(thread_pool[i], NULL);
    }
    pthread_attr_destroy(&attr);
}

type_t *allocate_real_vector(int l, int u)
{ /* Allocates a real vector of range [l..u]. */
    // void system-error (char *) ;
    type_t *p;
    p = (type_t *)malloc((unsigned)(u - l + 1) * sizeof(type_t));
    // if ( !p) system-error ("Failure in allocate-realgector 0 . "1 ;
    return p - 1;
}

type_t **allocate_real_matrix(int lr, int ur, int lc, int uc)
{
    /* Allocates a real matrix of range [lr..url [lc..ucl. */
    // void system_error(char *);
    int i;
    type_t **p;
    p = (type_t **)malloc((unsigned)(ur - lr + 1) * sizeof(type_t *));
    // if ( !p) system-error ("Failure in allocate-real-matrix 0 . ") ;
    p -= lr;
    for (i = lr; i <= ur; i++)
    {
        p[i] = (type_t *)malloc((unsigned)(uc - lc + 1) * sizeof(type_t));
        // if ( !p [i] ) system-error ("Failure in allocate-real-matrix 0 . ") ;
        p[i] -= lc;
    }
    return p;
}
void free_real_vector(type_t *v, int l)
{
    // Frees a real vector of range [l...u]
    free((void *)(v + l));
}

void free_real_matrix(type_t **m, int lr, int ur, int lc)
{
    /* Frees a real matrix of range [lr. .url [lc. .uc] . */
    int i;
    for (i = ur; i >= lr; i--)
    {
        free((void *)(m[i] + lc));
    }
    free((char *)(m + lr));
}