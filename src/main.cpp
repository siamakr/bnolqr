#include <Arduino.h>
#include <Wire.h>
#include <Adafruit_Sensor.h>
#include <Adafruit_BNO055.h>
#include <utility/imumaths.h>
#include <Adafruit_SPIDevice.h>
#include <Wire.h>
#include <Servo.h>
#include <Math.h>
#include <BasicLinearAlgebra.h>

#define PI 3.141592653589793238462643383279
#define MAX_ANGLE 120.0 //degrees
#define I_LIMIT 5    //degrees
#define INTERTIA 0.496753
#define COM_TO_TVC 0.425
#define MASS 2.372                    //Kg
#define MAX_ANGLE_SERVO 15           //Deg
#define SERVOMAX 35                   //Deg
#define SERVOMAXPWM 2300              //uSec
#define SERVOMIN -35                  //Deg
#define SERVOMINPWM 999               //uSec
#define TVC_MAX_ANGLE 24.5
#define Y_TORQUE_MAX 6.00
#define MAX_ROT_PER_STEP 2.16667/.01      //deg

//PITCH // max: 1862 //min: 1162       //mid:1512       //mid:1512
//ROLL: max:L 1725  //min:1207  //mid: 1412

#define PITCH_TVC_CENTER_PWM 1462     //uSec 
#define PITCH_TVC_MAX_PWM 1950        //uSec
#define PITCH_TVC_MIN_PWM 1020        //uSec
#define ROLL_TVC_CENTER_PWM 1500       //uSec
#define ROLL_TVC_MAX_PWM 1885          //uSec
#define ROLL_TVC_MIN_PWM 1060          //uSec
#define EDF_OFF_PWM 1080              //uSec
#define EDF_MIN_PWM 1100              //uSec
#define EDF_MAX_PWM 1800              //uSec
#define EDF_MAX_SUSTAINED_PWM 1730    //uSec
#define EDF_IDLE_PWM 1200             //uSec

#define EDF_THROTTLE_OFF_PWM 1000    
#define EDF_THROTTLE_MIN_PWM 1100     //uSec
#define EDF_THROTTLE_IDLE_PWM 1150    //.01 Kg Thrust
#define EDF_THROTTLE_SAFE_MAX_PWM 1480     //~1.1Kg thrust
#define EDF_THROTTLE_MAX_PWM 2000       //2KG thrust

#define Y_POLYVAL 0.0610
#define X_POLYVAL 0.0617

//VALUES SET FOR +-15º GIMBAL ANGLE FOR BOTH X AND Y
#define SERVO_INC_STEP_SIZE 10
#define SERVO_X_CENTER_US 1440
#define SERVO_Y_CENTER_US 1535
#define SERVO_X_MIN_US 1224
#define SERVO_Y_MIN_US 1160
#define SERVO_X_MAX_US 1892
#define SERVO_Y_MAX_US 1836

//SERVO-ANGLE TO TVC ANGLE POLYNOMIAL TRANSOFORMATOR 
#define X_P1 -29.2198405328854 
#define X_P2 1453.88991021228
#define Y_P1 -20.4054981741083 
#define Y_P2 1530.81204806643

//EDF MOTOR RPM TO THRUST TO PWM TRANSFORMATORS
#define RAD2N_P1 0.018566536813619    //Newtons to radians/s 
#define RAD2N_P2 -22.506778362213     //Newtons to Radians/s
#define RAD2PWM_P1 0.2719786528784    //Radians/s to PWM(us)  
#define RAD2PWM_P2 1010.29617153703   //Radians/s to PWM(us)
#define RPM_TO_OMEGA (2*PI/60)        //RPM to Radians/s
#define OMEGA_TO_RPM (60/(2*PI))      //Radians/s to RPM
#define GRAMS_TO_NEWTONS (9.8 / 1000) //Grams(g) to Newtons(N)

//MASS-MOMENT-OF-INERTIA OF VEHICLE
#define V_JXX 6.85563961956689
#define V_JYY 7.08910783725713
#define V_JZZ 0.0120276855157049
//MASS-MOMENT-OF-INERTIA OF EDF-PROP/MOTOR
#define EDF_JZZ 0.000174454522116462
//MASS-MOMENT-OF-INERTIA OF REACTION WHEEL 
#define RW_JZZ 0.00174245619164574

using namespace BLA;


int c = 0;
int count = 0;
float alpha = .04;      //Complimentary Filter Param
const int MPU = 0x68; // MPU6050 I2C address
float AccX, AccY, AccZ, GyroX, GyroY, GyroZ;
float pitch_buffer[3], roll_buffer[3], yaw_buffer[3];
float accAngleX, accAngleY, accAngleZ, gyroAngleX, gyroAngleY, gyroAngleZ;
float roll, pitch, yaw, rollK, pitchK, yawK, rollKal[2];

float curr_time, prev_time, dt;
float previousTime, currentTime, elapsedTime;
float GyroError[2], AccError[2];
float magX, magY, magZ;

Matrix<3,1> U; // Output vector
Matrix<6,1> error; // State error vector
Matrix<3,6> K; 
Matrix<6,1> REF; 
Matrix<6,1> X;

Matrix<4,1> U_red; // Output vector
Matrix<6,1> error_red; // State error vector
Matrix<4,6> K_red; 
Matrix<6,1> REF_red; 
Matrix<6,1> X_red;

Matrix<6,1> prev;

struct previous
{
  float gx; 
  float gy;
  float gz; 
  float u1; 
  float u2; 
  float u3;
};




float limit(float value, float min, float max);
float deg_to_rad(float deg);
float rad_to_deg(float rad);
void control_attitude(float r, float p, float y, float gx, float gy, float gz);
void control_attitude_red(float roll, float pitch, float yaw, float gx, float gy, float gz);
//void update_IMU(void);
//void IMU_init(void);
void calculate_IMU_error(void);
float LPF( float new_sample, float old_sample, float alpha );
float averaging(float new_sample, float old_sample);

/* This driver uses the Adafruit unified sensor library (Adafruit_Sensor),
   which provides a common 'type' for sensor data and some helper functions.
   To use this driver you will also need to download the Adafruit_Sensor
   library and include it in your libraries folder.
   You should also assign a unique ID to this sensor for use with
   the Adafruit Sensor API so that you can identify this particular
   sensor in any data logs, etc.  To assign a unique ID, simply
   provide an appropriate value in the constructor below (12345
   is used by default in this example).
   Connections
   ===========
   Connect SCL to analog 5
   Connect SDA to analog 4
   Connect VDD to 3.3-5V DC
   Connect GROUND to common ground
   History
   =======
   2015/MAR/03  - First release (KTOWN)
*/

/* Set the delay between fresh samples */
uint16_t BNO055_SAMPLERATE_DELAY_MS = 100;

// Check I2C device address and correct line below (by default address is 0x29 or 0x28)
//                                   id, address
Adafruit_BNO055 bno = Adafruit_BNO055(55, 0x28);

void printEvent(sensors_event_t* event);

void setup(void)
{
  Serial.begin(9600);
  Serial.println("Orientation Sensor Test"); Serial.println("");

  /* Initialise the sensor */
  if (!bno.begin())
  {
    /* There was a problem detecting the BNO055 ... check your connections */
    Serial.print("Ooops, no BNO055 detected ... Check your wiring or I2C ADDR!");
    while (1);
  }
  delay(1000);

U = {0,0,0};
error = {0,0,0,0,0,0};
REF = {0,0,0,0,0,0};
X = {0,0,0,0,0,0};

U_red = {0,0,0,0};
error_red = {0,0,0,0,0,0};
REF_red = {0,0,0,0,0,0};
X_red = {0,0,0,0,0,0};
prev = {0,0,0,0,0,0};


//K= {8.729696e-02,  0,  0,  6.313498e-01,  0,  0,  0,  8.729696e-02,  7.697951e-18,  0,  6.413067e-01,  -6.259864e-17,  0,  -1.356320e-17,  7.071068e-02,  0,  -1.730818e-17,  4.473222e+00};
K= {4.245702e-01,0,0,1.964807e-02,0,0,0,4.279157e-01,2.473468e-16,0,1.964807e-02,1.507959e-18,0,4.532035e-18,2.196117e+00,0,9.058981e-19,1.964807e-02};
//K_red = {-1.000000e+01,  0,  0,  0,  0,  -1.000000e+01,  0,  0,  0,  0,  3.684229e+03,  -1.796386e+00,  -6.604855e+00,  0,  0,  0,  0,  -6.713842e+00,  0,  0,  0,  0,  3.684355e+03,  -1.808265e+00,  0,  0,  5.509038e+03,  4.805412e+00,  0,  0,  3.006950e+03,  2.654509e+00};
}

void loop(void)
{
  //could add VECTOR_ACCELEROMETER, VECTOR_MAGNETOMETER,VECTOR_GRAVITY...

     imu::Vector<3> euler = bno.getVector(Adafruit_BNO055::VECTOR_EULER);
     imu::Vector<3> gyro = bno.getVector(Adafruit_BNO055::VECTOR_GYROSCOPE);

  /* Display the floating point data */
  Serial.print("X: ");
  Serial.print(euler.x());
  Serial.print("\t\t");
  Serial.print(" Y: ");
  Serial.print(euler.y());
  Serial.print("\t\t");
  Serial.print(" Z: ");
  Serial.print(euler.z());
  Serial.print("\t\t");
  Serial.println("");
   //control_attitude(roll, pitch, yaw, GyroX, GyroY, GyroZ);
   control_attitude(deg_to_rad(euler.z()), deg_to_rad(euler.y()), deg_to_rad(euler.x()), deg_to_rad(gyro.z()), deg_to_rad(gyro.y()), deg_to_rad(gyro.x()));
  // control_attitude(deg_to_rad(deg_to_rad(euler.z())), deg_to_rad(euler.y()), deg_to_rad(euler.x()), deg_to_rad(gyro.z()), deg_to_rad(gyro.y()), deg_to_rad(gyro.x()));
  //printEvent(&orientationData);
  //printEvent(&angVelocityData);
  //printEvent(&linearAccelData);
  //printEvent(&magnetometerData);
 // printEvent(&accelerometerData);
  //printEvent(&gravityData);


  //delay(BNO055_SAMPLERATE_DELAY_MS);
}
void updatebno(void){
    //could add VECTOR_ACCELEROMETER, VECTOR_MAGNETOMETER,VECTOR_GRAVITY...
  sensors_event_t orientationData , angVelocityData , linearAccelData, magnetometerData, accelerometerData, gravityData;
  bno.getEvent(&orientationData, Adafruit_BNO055::VECTOR_EULER);
  bno.getEvent(&angVelocityData, Adafruit_BNO055::VECTOR_GYROSCOPE);
  bno.getEvent(&linearAccelData, Adafruit_BNO055::VECTOR_LINEARACCEL);
  bno.getEvent(&magnetometerData, Adafruit_BNO055::VECTOR_MAGNETOMETER);
  bno.getEvent(&accelerometerData, Adafruit_BNO055::VECTOR_ACCELEROMETER);
  bno.getEvent(&gravityData, Adafruit_BNO055::VECTOR_GRAVITY);



}

void calibratebno(void){
  
  uint8_t system, gyro, accel, mag = 0;
  bno.getCalibration(&system, &gyro, &accel, &mag);
  //Serial.println();
  //Serial.print("Calibration: Sys=");
  //Serial.print(system);
  // Serial.print(" Gyro=");
  // Serial.print(gyro);
  // Serial.print(" Accel=");
  // Serial.print(accel);
  // Serial.print(" Mag=");
  // Serial.println(mag);

}
void control_attitude(float r, float p, float y, float gx, float gy, float gz){
  //X < roll, pitch, yaw, gx, gy, gz; 
  X(1) = r;
  X(2) = p;
  X(3) = 0;
  X(4) = LPF(gx,prev(1), .5);
  X(5) = LPF(gy,prev(2), .5);
  X(6) = LPF(gz, prev(3), 0.3);
  
  error = X - REF; 
  U = -K * error; 

// float u1temp = averaging(U(1), prev(4));
// float u2temp = averaging(U(2), prev(5));
float u1temp = rad_to_deg(LPF(U(1), prev(4), .3));
float u2temp = rad_to_deg(LPF(U(2), prev(5), 0.3)); 



 Serial.print("u1: ");
 Serial.print(U(1));
 Serial.print("\t");
 Serial.print("u2: ");
 Serial.print(U(2));
 Serial.print("\t");
 Serial.print("u3: ");
 Serial.print(U(3));
 Serial.print("\t");
 Serial.print("servo1: ");
 Serial.print(u1temp);
 Serial.print("\t");
 Serial.print("servo2: ");
 Serial.print(u2temp);
 Serial.print("\t");
//  Serial.println("");

//shift arrays 
prev(1) = X(1) ; 
prev(2) = X(2) ; 
prev(3) = X(3) ; 
prev(4) = LPF(U(1), prev(4), .3) ; 
prev(5) = LPF(U(2), prev(5), 0.3) ; 
prev(6) = 0; 


}
void control_attitude_red(float roll, float pitch, float yaw, float gx, float gy, float gz){
  //X < roll, pitch, yaw, gx, gy, gz; 
  X_red(1) = roll;
  X_red(2) = pitch;
  X_red(3) = 0;
  X_red(4) = gx;
  X_red(5) = gy;
  X_red(6) = gz;
  
  error_red = X_red - REF_red; 
  U_red = -K_red * error_red; 

//  U_red(1) = limit(U_red(1), -15, 15);
//  U_red(2) = limit(U_red(2), -15, 15); 

 Serial.print("u1: ");
 Serial.print(U_red(3));
 Serial.print("\t");
  Serial.print("edf2: ");
 Serial.print(U_red(4));
 Serial.print("\t");
 Serial.print("servo1: ");
 Serial.print(U_red(1));
 Serial.print("\t");
 Serial.print("servo2: ");
 Serial.print(U_red(2));
 Serial.print("\t");
//  Serial.println("");

}

float rad_to_deg(float rad){

  return rad * 180 / PI;
}

float deg_to_rad(float deg){

  return deg * PI / 180;
}

float limit(float value, float min, float max){

  if(value > max ) value = max; 
  if(value < min ) value = min; 

  return value; 
}
void printEvent(sensors_event_t* event) {
  double x = -1000000, y = -1000000 , z = -1000000; //dumb values, easy to spot problem
  if (event->type == SENSOR_TYPE_ACCELEROMETER) {
    Serial.print("Accl:");
    x = event->acceleration.x;
    y = event->acceleration.y;
    z = event->acceleration.z;
  }
  else if (event->type == SENSOR_TYPE_ORIENTATION) {
    Serial.print("Orient:");
    x = event->orientation.x;
    y = event->orientation.y;
    z = event->orientation.z;
  }
  else if (event->type == SENSOR_TYPE_MAGNETIC_FIELD) {
    Serial.print("Mag:");
    x = event->magnetic.x;
    y = event->magnetic.y;
    z = event->magnetic.z;
  }
  else if (event->type == SENSOR_TYPE_GYROSCOPE) {
    Serial.print("Gyro:");
    x = event->gyro.x;
    y = event->gyro.y;
    z = event->gyro.z;
  }
  else if (event->type == SENSOR_TYPE_ROTATION_VECTOR) {
    Serial.print("Rot:");
    x = event->gyro.x;
    y = event->gyro.y;
    z = event->gyro.z;
  }
  else if (event->type == SENSOR_TYPE_LINEAR_ACCELERATION) {
    Serial.print("Linear:");
    x = event->acceleration.x;
    y = event->acceleration.y;
    z = event->acceleration.z;
  }
  else if (event->type == SENSOR_TYPE_GRAVITY) {
    Serial.print("Gravity:");
    x = event->acceleration.x;
    y = event->acceleration.y;
    z = event->acceleration.z;
  }
  else {
    Serial.print("Unk:");
  }

  Serial.print("  \tx= ");
  Serial.print(x);
  Serial.print("  \ty= ");
  Serial.print(y);
  Serial.print("  \tz= ");
  Serial.println(z);
}

float LPF( float new_sample, float old_sample, float alpha ){
    return ((alpha * new_sample) + (1.0-alpha) * old_sample);  
}

float averaging(float new_sample, float old_sample){
  return ((new_sample + old_sample)/2);
}
