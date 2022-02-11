# fpga-6dof-inverse-kinematics

*This implementation will be on Xilinx Artix-7 100T FPGA (Digilent Arty A7-100T Board) with Vivado using Verilog and host the RISC-V processor.* 

It introduces a hardware synthesized inverse kinematics algorithm that controls 6 degrees of freedom robotic arms and a sensor-based automatic control system. All systems are hardware-synthesized and work simultaneously. (2022-06-25~ (If I have some free time, I'll start working on it.)) 

Before starting the project, I implement an algorithm that seeks numerical solutions as long as I have free time. I will port this shape to Verilog. To facilitate Verilog porting, I'm writing as hardware-friendly code as possible, focusing on the C language. ROS is not used in this project because it aims to control in a low-power embedded environment.

# Mechanical Foundation

## Relationship between links.
![graphical explanation](https://user-images.githubusercontent.com/71680670/153541844-dac2ad61-db49-498f-80b7-cd1af2921975.PNG)

## Geometric characteristics of 6dof inverse kinematics
![explanation of codegen_m](https://user-images.githubusercontent.com/71680670/153542220-a1177275-32f7-4af8-a11f-ffe19a3725be.PNG)

## Multivariable Newton's method for solving nonlinear kinematics problems
![how to solve the inverse kinematics problem](https://user-images.githubusercontent.com/71680670/153542398-7351df60-5a76-47dd-9ec4-42cf083dffab.PNG)
