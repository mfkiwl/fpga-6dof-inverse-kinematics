#include "pwm_servo_control.h"

uint16_t pulse_width(long double angle)
{
  uint16_t analog_value;
  long double pulse_wide = MIN_PULSE_WIDTH + ((angle - MIN_ANLGE) / (MAX_ANGLE - MIN_ANLGE)) * (MAX_PULSE_WIDTH - MIN_PULSE_WIDTH);
  analog_value = (uint16_t)(round((pulse_wide / 1000000.0) * FREQUENCY * 4096));
  return analog_value;
}

void rot_servo(Adafruit_PWMServoDriver *pwm, int servo_num, angle_disp_t disp, deg_per_sec_t default_anglular_velocity)
{
  if (disp.th_i != UNKNOWN_ANGLE)
  {
    double tot_delay_time = disp.taken_time - abs(disp.th_f - disp.th_i) / default_anglular_velocity;
    if (tot_delay_time > 0)
    {
      angle_t dth = (disp.th_f - disp.th_i) / SPEED_ADJUSTMENT_RESOLUTION;
      angle_t cur_th;
      time_t delay_millis = round((tot_delay_time * 1000.0) / SPEED_ADJUSTMENT_RESOLUTION);
      for (cur_th = disp.th_i; dth * (disp.th_f-cur_th) > 0; cur_th += dth)
      {
        pwm->setPWM(servo_num, 0, pulse_width(cur_th));
        delay(delay_millis);
      }
      // end of rotation and then return to main 
      exit(0);
    }
  }
  pwm->setPWM(servo_num, 0, pulse_width(disp.th_f));
}

void rot_mult_servo(Adafruit_PWMServoDriver *pwm, int servo_num_arr, angle_disp_t disp_arr, int arr_len)
{
  //it will be implemented by POSIX thread  
}