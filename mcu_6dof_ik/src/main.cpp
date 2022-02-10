#include <pwm_servo_control.h>
#include <inverse_kinematics.h>
#include <stdlib.h>
#include <stdio.h>

Adafruit_PWMServoDriver pwm = Adafruit_PWMServoDriver(I2C_DEFAULT_ADDR);
mat_element_t init_vector[6] = {0.0};
mat_element_t fin_vector[6] = {0.0};

void setup()
{
  Serial.begin(9600);
  Serial.println("6 dof rotbotic arm inverse kinematics");
  pwm.begin();
  pwm.setPWMFreq(FREQUENCY);
  init_servo_location(&pwm, init_vector);
}

bool is_first = true;

void memcpy(mat_element_t* a, mat_element_t* b, size_t c)
{
  int len = c / sizeof(mat_element_t);
  for (int i = 0; i < len; i++)
  {
    a[i] = b[i];
  }
}
void loop()
{
  mat_element_t vec0[6];
  mat_element_t vec1[6] = {90, 90, 90, 90, 90, 90};
  mat_element_t vec2[6] = {180, 135, 0, 90, 0, 90};

  init_servo_location(&pwm, init_vector);
  memcpy(vec0, init_vector, sizeof(mat_element_t) * 6);
  memcpy(fin_vector, vec1, sizeof(mat_element_t) * 6);
  rot_mult_servo(&pwm, init_vector, fin_vector);

  memcpy(init_vector, vec1, sizeof(mat_element_t) * 6);
  memcpy(fin_vector, vec2, sizeof(mat_element_t) * 6);
  rot_mult_servo(&pwm, init_vector, fin_vector);

  memcpy(init_vector, vec2, sizeof(mat_element_t) * 6);
  memcpy(fin_vector, vec0, sizeof(mat_element_t) * 6);
  rot_mult_servo(&pwm, init_vector, fin_vector);

  delay(1000);
  exit(0);
}