#ifndef PWM_SERVO_CONTROL_H
#define PWM_SERVO_CONTROL_H

#include <Arduino.h>
#include <time.h>
#include <Wire.h>
#include "..\..\.pio\libdeps\esp-wrover-kit/Adafruit PWM Servo Driver Library\Adafruit_PWMServoDriver.h"
#include <math.h>
#include "../inverse_kinematics/inverse_kinematics.h"
#include <stdbool.h>
#include <stdlib.h>

#define MIN_PULSE_WIDTH 500 
#define MAX_PULSE_WIDTH 2500 
#define DEFAULT_PULSE_WIDTH 1500
#define MIN_ANLGE 0 
#define MAX_ANGLE 180 
#define FREQUENCY 50
#define I2C_DEFAULT_ADDR 0x40
#define SPEED_ADJUSTMENT_RESOLUTION 1000
#define UNKNOWN_ANGLE -1

typedef double deg_per_sec_t;
typedef struct _angle_disp_s
{
    angle_t th_i;
    angle_t th_f;
    time_t taken_time;
} angle_disp_t;

uint16_t pulseWidth(long double angle);
void rot_servo(Adafruit_PWMServoDriver* pwm, int servo_num, angle_disp_t dth, deg_per_sec_t default_vel);
void rot_mult_servo(Adafruit_PWMServoDriver *pwm, int servo_num_arr, angle_disp_t disp_arr, int arr_len);

#endif
