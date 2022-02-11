#ifndef LINALG_H
#define LINALG_H

//#include "inverse_kinematics.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <stdbool.h>


#define MAX_LEN 6
#define NUM_T 1


#define l1 1
#define l2 1
#define l3 1
#define l4 1
#define l5 1
#define l6 1

#define DONE 0
#define DIV_BY_ZERO -1

typedef double type_t;

typedef struct arguments
{
    int start_row;
    int end_row;
} arguments;


void multiply(int N,  type_t mat_a[][MAX_LEN], type_t mat_b[][MAX_LEN], type_t mat_c[][MAX_LEN]);
void get_jacobian_with_pthread(type_t jacobian[][6], type_t cur_angle[6]);
void get_jacobian(type_t jacobian[][6], type_t cur_angle[6]);
int inverse_matrix(int dim, type_t inv_a[][MAX_LEN], type_t a[][MAX_LEN]);




#endif
