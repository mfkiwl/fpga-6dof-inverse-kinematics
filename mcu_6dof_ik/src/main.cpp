#include <pwm_servo_control.h>
#include <inverse_kinematics.h>
#include <stdlib.h>
#include <stdio.h>

Adafruit_PWMServoDriver pwm = Adafruit_PWMServoDriver(I2C_DEFAULT_ADDR);

void setup()
{
  pwm.begin();
  pwm.setPWMFreq(FREQUENCY);
}

void loop()
{
  // for testing after the base method operation in the lib directory.
}