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
#define COM_TO_TVC 0.1335
#define MASS 2.458                    //Kg
#define MAX_ANGLE_SERVO 15           //Deg
#define X 0
#define Y 1
#define Z 2
#define SERVO_MAX_SEC_PER_DEG 0.001667
#define DT .01
#define SERVO_MAX_DEGREE_PER_DT 12

//PITCH // max: 1862 //min: 1162       //mid:1512       //mid:1512
//ROLL: max:L 1725  //min:1207  //mid: 1412

#define PITCH_TVC_CENTER_PWM 1462     //uSec 
#define PITCH_TVC_MAX_PWM 1950        //uSec
#define PITCH_TVC_MIN_PWM 1020        //uSec
#define ROLL_TVC_CENTER_PWM 1500       //uSec
#define ROLL_TVC_MAX_PWM 1885          //uSec
#define ROLL_TVC_MIN_PWM 1060          //uSec
#define EDF_OFF_PWM 1000              //uSec
#define EDF_MIN_PWM 1500              //uSec
#define EDF_MAX_PWM 2000              //uSec
#define EDF_MAX_SUSTAINED_PWM 1730    //uSec
#define EDF_IDLE_PWM 1600             //uSec

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
#define YY_P1 -.204054981741083 
#define YY_P2 1530.81204806643

//EDF MOTOR RPM TO THRUST TO PWM TRANSFORMATORS
#define RAD2N_P1 0.018566536813619    //Newtons to radians/s 
#define RAD2N_P2 -22.506778362213     //Newtons to Radians/s
#define RAD2PWM_P1 0.2719786528784    //Radians/s to PWM(us)  
#define RAD2PWM_P2 1010.29617153703   //Radians/s to PWM(us)
#define RPM_TO_OMEGA (2*PI/60)        //RPM to Radians/s
#define OMEGA_TO_RPM (60/(2*PI))      //Radians/s to RPM
#define GRAMS_TO_NEWTONS (9.8 / 1000) //Grams(g) to Newtons(N)



//MASS-MOMENT-OF-INERTIA OF VEHICLE
#define V_JXX 0.0058595
#define V_JYY 0.0058595
#define V_JZZ 0.0120276855157049
//MASS-MOMENT-OF-INERTIA OF EDF-PROP/MOTOR
#define EDF_JZZ 0.000174454522116462
//MASS-MOMENT-OF-INERTIA OF REACTION WHEEL 
#define RW_JZZ 0.00174245619164574

using namespace BLA;


// int c = 0;
// int count = 0;
// float alpha = .04;      //Complimentary Filter Param
// const int MPU = 0x68; // MPU6050 I2C address
// float AccX, AccY, AccZ, GyroX, GyroY, GyroZ;
// float pitch_buffer[3], roll_buffer[3], yaw_buffer[3];
// float accAngleX, accAngleY, accAngleZ, gyroAngleX, gyroAngleY, gyroAngleZ;
// float roll, pitch, yaw, rollK, pitchK, yawK, rollKal[2];

float curr_time, prev_time, dt;
float previousTime, currentTime, elapsedTime;
float sensor_timer, control_timer; 
//float GyroError[2], AccError[2];


//Matrix<3,1> U = {0,0,0}; // Output vector
//Matrix<6,1> error = {0,0,0,0,0,0}; // State error vector
Matrix<3,6> K =   {.9119,    0.0000,   -0.0000,   -0.52107,    0.0000,   -0.0000,
                  0.0000,    0.72119,   -0.0000,   -0.0000,    -0.50107,   -0.0000,
                  0.0000,   -0.0000,    24.180,   -0.0000,    0.0000,   0.471857};
Matrix<6,1> REF = {0, 0, 0, 0, 0, 0}; 
//Matrix<6,1> Xs = {0, 0, 0, 0, 0, 0};

Matrix<3,1> U_red = {0,0,0}; // Output vector
Matrix<6,1> error_red = {0,0,0,0,0,0}; // State error vector
Matrix<3,6> K_red = {0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0}; 
Matrix<6,1> REF_red = {0, 0, 0, 0, 0, 0}; 
Matrix<6,1> Xs_red = {0, 0, 0, 0, 0, 0};

Matrix<6,1> prev;

Servo sx; 
Servo sy; 
Servo edf;

typedef struct
{
  float roll, pitch, yaw; 
  float gx, gy, gz;
  float ax, ay, az; 
  float lax, lay, laz; 
  float vx, vy, vz; 
  float grx, gry, grz; 
  float x, y, z; 
  float magx, magy, magz; 
  float u1, u2, u3, u4; 
  imu::Vector<3> trq; 
  imu::Vector<3> T;
} curr_t;

typedef struct
{
  float roll, pitch, yaw; 
  float gx, gy, gz;
  float ax, ay, az; 
  float lax, lay, laz; 
  float vx, vy, vz; 
  float grx, gry, grz; 
  float x, y, z; 
  float magx, magy, magz; 
  float u1, u2, u3, u4; 
  imu::Vector<3> trq; 
  imu::Vector<3> T;
} prev_t;

typedef struct
{
  float roll, pitch, yaw; 
  float gx, gy, gz;
  float ax, ay, az; 
  float lax, lay, laz; 
  float vx, vy, vz; 
  float grx, gry, grz; 
  float x, y, z; 
  float magx, magy, magz; 
  float u1, u2, u3, u4;
  imu::Vector<3> trq; 
  imu::Vector<3> T; 
} pprev_t;

float T[2], torque[2];
  float servoang1, Tmag;
  float servoang2 ;



// function definitions (will need to add these to the header file but this is just for initial testing puroses to keep everything in one file)
float limit(float value, float min, float max);
float d2r(float deg);
float r2d(float rad);
void control_attitude(float r, float p, float y, float gx, float gy, float gz);
void control_attitude_red(float roll, float pitch, float yaw, float gx, float gy, float gz);
void samplebno(void);
//void update_IMU(void);
//void IMU_init(void);
void calculate_IMU_error(void);
float LPF( float new_sample, float old_sample, float alpha );
float averaging(float new_sample, float old_sample);
void writeXservo(float angle);
void writeYservo(float angle);
void init_servos(void);
void writeEDF(float Ft);
float ft2omega(float Ft);
int omega2pwm(float omega);
void emergency_check(float r, float p);
void init_edf(void);
float servo_rate_check(float new_sample, float old_sample);


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
uint16_t BNO055_SAMPLERATE_DELAY_MS = 10;

// Check I2C device address and correct line below (by default address is 0x29 or 0x28)
//                                   id, address
Adafruit_BNO055 bno = Adafruit_BNO055(55, 0x28);

void printEvent(sensors_event_t* event);
//states current, t-1 and t-2 (for filtering)
curr_t c;
prev_t p; 
pprev_t pp; 

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
  init_servos();


 K= {.9119,    0.0000,   -0.0000,   -0.52107,    0.0000,   -0.0000,
      0.0000,    0.72119,   -0.0000,   -0.0000,    -0.50107,   -0.0000,
      0.0000,   -0.0000,    24.180,   -0.0000,    0.0000,   0.471857};

init_servos();
init_edf();

}

void loop(void)
{
  //could add VECTOR_ACCELEROMETER, VECTOR_MAGNETOMETER,VECTOR_GRAVITY...
  if(micros() - sensor_timer >= 20000){
    sensor_timer = micros();
    samplebno(); 
  }
  
  if(micros() - control_timer >= 5000){
    control_timer = micros();
    control_attitude(d2r(c.roll), d2r(c.pitch), d2r(c.yaw), d2r(c.gx), d2r(c.gy), d2r(c.gz)); 
  }



   control_attitude(d2r(euler.z()), d2r(euler.y()), d2r(euler.x()), d2r(gyro.z()), d2r(gyro.y()), d2r(gyro.x()));

  emergency_check(euler.z(), euler.y());


delay(BNO055_SAMPLERATE_DELAY_MS);
}
void updatebno(void){

  imu::Vector<3> euler = bno.getVector(Adafruit_BNO055::VECTOR_EULER);
  imu::Vector<3> gyro = bno.getVector(Adafruit_BNO055::VECTOR_GYROSCOPE);
  imu::Vector<3> accel = bno.getVector(Adafruit_BNO055::VECTOR_ACCELEROMETER);
  imu::Vector<3> linaccel = bno.getVector(Adafruit_BNO055::VECTOR_LINEARACCEL);
  imu::Vector<3> grav = bno.getVector(Adafruit_BNO055::VECTOR_GRAVITY);
  imu::Vector<3> mag = bno.getVector(Adafruit_BNO055::VECTOR_MAGNETOMETER);

//Shift the time state vectors from t-1 to t-2, from t-0 to t-1; 

//shift from "current" (c) state to "previous" (p) state 
  p.roll = c.roll; p.pitch = c.pitch; p.yaw =  c.yaw;
  p.gx =  c.gx; p.gy =  c.gy; p.gz =  c.gz;
  p.ax =  c.ax; p.ay =  c.ay; p.az =  c.az;
  p.lax =  c.lax; p.lay =  c.lay; p.laz =  c.laz;
  p.grx =  c.grx; p.gry = c.gry; p.grz =  c.grz;
  p.magx =  c.magx; p.magy =  c.magy; p.magz =  c.magz;

 //shift from "previous state " (p) to "previous-previous state" (pp) 
  pp.roll = p.roll; pp.pitch = p.pitch; pp.yaw =  p.yaw; 
  pp.gx =  p.gx; pp.gy =  p.gy; pp.gz =  p.gz;
  pp.ax =  p.ax; pp.ay =  p.ay; pp.az =  p.az;
  pp.lax =  p.lax; pp.lay =  p.lay; pp.laz =  p.laz;
  pp.grx =  p.grx; pp.gry =  p.gry; pp.grz =  p.grz;
  pp.magx =  p.magx; pp.magy =  p.magy; pp.magz =  p.magz;
  
  //Sample IMU and put values into "current state" (c) struct
  c.roll = euler.z(); c.pitch = euler.y(); c.yaw = euler.x();
  c.gx = gyro.z(); c.gy = gyro.y(); c.gz = gyro.x();
  c.ax = accel.z(); c.ay = accel.y(); c.az = accel.x();
  c.lax = linaccel.z(); c.lay = linaccel.y(); c.laz = linaccel.x();
  c.grx = grav.z(); c.gry = grav.y(); c.grz = grav.x();
  c.magx = mag.z(); c.magy = mag.y(); c.magz = mag.x();



  


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
void control_attitude(float Roll, float Pitch, float Yaw, float GX, float GY, float GZ){

Matrix<3,1> U = {0,0,0};              // Output vector
Matrix<6,1> error = {0,0,0,0,0,0};    // State error vector
Matrix<6,1> Xs = {0, 0, 0, 0, 0, 0};  // Temporary state vector
  
  //shifting the control outputs to previous time steps
  //shift from T-1 (p) to T-2 (pp)
  pp.u1 = p.u1; pp.u2 = p.u2; pp.u3 = p.u3; pp.u4 = p.u4;
  pp.trq.x() = p.trq.x(); pp.trq.y() = p.trq.y(); pp.trq.z() = p.trq.z();
  pp.T.x() = p.T.x(); pp.T.y() = p.T.y(); pp.T.z() = p.T.z();
  //shift from T (c) to T-1 (p)
  p.u1 = c.u1; p.u2 = c.u2; p.u3 = c.u3; p.u4 = c.u4;
  p.trq.x() = c.trq.x(); p.trq.y() = c.trq.y(); p.trq.z() = c.trq.z();
  p.T.x() = c.T.x(); p.T.y() = c.T.y(); p.T.z() = c.T.z();
  
  

  //load Xs state vector with measured values 
  Xs(0) = Roll;
  Xs(1) = Pitch;
  Xs(2) = 1;
  Xs(3) = LPF(GX, prev(0), .94);
  Xs(4) = LPF(GY, prev(1), .94);
  Xs(5) = 0;
  //Xs = {r, p, 0, gx, gy, 0};
  error = Xs-REF; 
  U = -K * error; 

  //load desired torque vector
  c.trq.x() = U(2)*sin(U(0))*COM_TO_TVC;
  c.trq.y() = U(2)*sin(U(1))*COM_TO_TVC;
  c.trq.z() = U(2);

  //load new Thrust Vector from desired torque
  c.T.x() = trq.x()/COM_TO_TVC; 
  c.T.y() = -trq.y()/COM_TO_TVC; 
  c.T.z() = U(2);           //constant for now, should be coming from position controller 

  //calculate thrust magnitude 
  c.u3 = sqrt(pow(T.x(),2) + pow(T.y(),2) + pow(T.z(),2)); 
  
  //get servo angles from thrust vector 
  servoang1 = asin(T.x()/(Tmag - pow(T.y(),2)));
  servoang2 = asin(T.y()/Tmag);

// float u1temp = averaging(U(1), prev(4));
// float u2temp = averaging(U(2), prev(5));
c.u1 = limit(LPF(U(0), p.u1, .45), d2r(-15),d2r(15));
c.u2 = limit(LPF(U(1), p.u2, .45), d2r(-15),d2r(15)); 

 writeXservo(r2d(c.u1));
 writeYservo(r2d(c.u2));
 //writeEDF(Tmag);

//  Serial.print("r:  ");
//  Serial.print(r2d(r));
//  Serial.print("\t");
//  Serial.print("p:  ");
//  Serial.print(r2d(p));
//  Serial.print("\t");
//  Serial.print("gx:  ");
//  Serial.print(r2d(gx));
//  Serial.print("\t");
//  Serial.print("gy:  ");
//  Serial.print(r2d(gy));
//  Serial.print("\t");
//  Serial.print("ang1: ");
//  Serial.print(r2d(u1temp));
//  Serial.print("\t");
//  Serial.print("ang2: ");
//  Serial.print(r2d(u2temp));
//  Serial.print("\t");
//  Serial.print("Tmag: ");
//  Serial.print(Tmag);
//  Serial.print("\t");
//  Serial.println("\t\t");

//shift arrays 




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

float r2d(float rad){

  return rad * 180 / PI;
}

float d2r(float deg){

  return deg * PI / 180;
}

float limit(float value, float min, float max){

  if(value >= max ) value = max; 
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

void init_servos(void){
  sx.attach(11);
  sy.attach(10);

  //Zero Servos 
  writeXservo(0);
  writeYservo(0);
}

void init_edf(void){
  edf.attach(9);
  delay(20);

  //initialize the edf motor and run it at 1500us for 5 seconds to initialize the govenor mode to linearize the throttle curve. 
  edf.writeMicroseconds(EDF_OFF_PWM); 
  delay(200);

  //go to 1500 and wait 5 seconds
  edf.writeMicroseconds(EDF_MIN_PWM);
  delay(5000);
}

void writeXservo(float angle){
  //map angle in degrees to pwm value for servo 
  int pwm = round((1000 * angle) * (X_P1/1000) + X_P2); 
  sx.writeMicroseconds(pwm); 
//   Serial.print("pwmx: ");
//  Serial.print(pwm);
//  Serial.print("\t");
}

void writeYservo(float angle){
  //map angle in degrees to pwm value for servo 
  int pwm = round((1000 * angle) * (Y_P1/1000) + Y_P2); // using polynomial regression coefficients to map tvc angle to pwm vals
  sy.writeMicroseconds(pwm); 
//   Serial.print("pwmy: ");
//  Serial.print(pwm);
//  Serial.print("\t");
}

void writeEDF(float Ft){
  float omega = (Ft - RAD2N_P2)/RAD2N_P1; 
  int pwm =round(omega * RAD2PWM_P1 + RAD2PWM_P2); 
//   Serial.print("pwmedf: ");
//  Serial.print(pwm);
//  Serial.print("\t");
  //write to edf servo  
}

float ft2omega(float Ft){
  return (Ft - RAD2N_P2)/RAD2N_P2;
}

int omega2pwm(float omega){
  return (int) round(omega * RAD2PWM_P1 + RAD2PWM_P2);
}

void emergency_check(float r, float p){
  if(r >= 45 || r <=-45 || p >= 45 || p <= -45){
    writeEDF(0); 
    writeXservo(0);
    writeYservo(0);
    while(1);
  }
}

float servo_rate_check(float new_sample, float old_sample){
  float temp; 
  if(new_sample > old_sample + d2r(SERVO_MAX_DEGREE_PER_DT)) temp = old_sample + d2r(SERVO_MAX_DEGREE_PER_DT);
  if(new_sample <= old_sample - d2r(SERVO_MAX_DEGREE_PER_DT)) temp = old_sample - d2r(SERVO_MAX_DEGREE_PER_DT);
  return temp; 

}

