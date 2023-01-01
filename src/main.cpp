#include <Arduino.h>
#include <Wire.h>
#include <Adafruit_Sensor.h>
#include <Adafruit_BNO055.h>
#include <utility/imumaths.h>
#include <Adafruit_SPIDevice.h>
#include <Servo.h>
#include <Math.h>
#include <BasicLinearAlgebra.h>
#include <stdint.h>
#include <LIDARLite_v3HP.h>
#include <float.h>
#include <String.h>

// Arduino Pin assignments
#define XSERVO_PIN 11
#define YSERVO_PIN 10
#define RW_PIN 6
#define EDF_PIN 9

#define COM_TO_TVC 0.1335
#define MASS 2.458                    //Kg
#define MAX_ANGLE_SERVO 15           //Deg
// #define X 0
// #define Y 1
// #define Z 2
//#define SERVO_MAX_SEC_PER_DEG 0.001667      //dt/.001667 |(dt=0.1) = 6º
#define SERVO_MAX_SEC_PER_DEG 0.003333f      //dt/.003333 |(dt=0.1) = 3º
#define SERVO_MAX_SEC_PER_RAD 0.0095496f      //dt/.003333 |(dt=0.1) = 3º
//#define DT .01
#define SERVO_MAX_DEGREE_PER_DT 12
#define MAX_VEHICLE_ANGLE_DEG 35.00f
#define DEADBAND_ANGLE_DEG 0.001f
#define SERVO_ANG_TO_TVC_ANG 3.00f

#define PITCH_TVC_CENTER_PWM 1462     //uSec 
#define PITCH_TVC_MAX_PWM 1950        //uSec
#define PITCH_TVC_MIN_PWM 1020        //uSec
#define ROLL_TVC_CENTER_PWM 1500       //uSec
#define ROLL_TVC_MAX_PWM 1885          //uSec
#define ROLL_TVC_MIN_PWM 1060          //uSec
#define EDF_OFF_PWM 900              //uSec
#define EDF_MIN_PWM 1500              //uSec
#define EDF_MAX_PWM 2000              //uSec
#define EDF_MAX_SUSTAINED_PWM 1730    //uSec
#define EDF_IDLE_PWM 1600             //uSec

//VALUES SET FOR +-15º GIMBAL ANGLE FOR BOTH X AND Y
#define SERVO_INC_STEP_SIZE 10
#define SERVO_X_CENTER_US 1440
#define SERVO_Y_CENTER_US 1535
#define SERVO_X_MIN_US 1224
#define SERVO_Y_MIN_US 1160
#define SERVO_X_MAX_US 1892
#define SERVO_Y_MAX_US 1836

//SERVO-ANGLE TO TVC ANGLE POLYNOMIAL TRANSOFORMATOR 
#define X_P1 -29.2198405328854f 
#define X_P2 1453.88991021228f
#define Y_P1 -20.4054981741083f 
#define Y_P2 1530.81204806643f
#define YY_P1 -.204054981741083f 
#define YY_P2 1530.81204806643f

//EDF MOTOR RPM TO THRUST TO PWM TRANSFORMATORS
#define RAD2N_P1 0.018566536813619f    //Newtons to radians/s 
#define RAD2N_P2 -22.506778362213f     //Newtons to Radians/s
#define RAD2PWM_P1 0.2719786528784f    //Radians/s to PWM(us)  
#define RAD2PWM_P2 1010.29617153703f   //Radians/s to PWM(us)
#define RPM_TO_OMEGA (2.0f*PI/60.0f)        //RPM to Radians/s
#define OMEGA_TO_RPM (60.0f/(2.0f*PI))      //Radians/s to RPM
#define GRAMS_TO_NEWTONS (9.80f / 1000.00f) //Grams(g) to Newtons(N)

//MASS-MOMENT-OF-INERTIA OF VEHICLE
#define V_JXX 0.0058595f
#define V_JYY 0.0058595f
#define V_JZZ 0.01202768f
//MASS-MOMENT-OF-INERTIA OF EDF-PROP/MOTOR
#define EDF_JZZ 0.0001744f
//MASS-MOMENT-OF-INERTIA OF REACTION WHEEL 
#define RW_JZZ 0.00174245f
#define DT_MSEC 5.00f
#define DT_SEC DT_MSEC/1000


using namespace BLA;

typedef struct
{
float roll, pitch, yaw; 
float gx, gy, gz;
float ax, ay, az;
float axw, ayw, azw;
float vx, vy, vzi, vzl, vz; 
float x, y, zb, z;
float u1, u2, u3, u4; 
} data_prev_t;

typedef struct 
{
  float roll, pitch, yaw;
  float gx, gy, gz;
  float ax, ay, az; 
  float axw, ayw, azw;
  float x, y, z, zraw; 
  float vx, vy, vzi, vzl, vz;
  float u1, u2, u3, u4, Tz, Tmag, ang1, ang2; 
}data_current_t;

typedef struct 
{
  float vx, vy, vz; 
  float x, y, z; 
}data_estimator_t;

/*
// Temporary IMU vectors, will be replaced with typdef structs
imu::Vector<3> euler ;
imu::Vector<3> gyro ;
imu::Vector<3> ef_r;
imu::Vector<3> gf_r;
imu::Vector<3> gflpr_r;
*/
/////////////////////// ATTITUDE ONLY START ///////////////////////
Matrix<3,1> U = {0,0,0}; // Output vector
Matrix<6,1> error {0,0,0,0,0,0}; // State error vector
Matrix<6,1> REF = {0,0,0,0,0,0}; 
Matrix<6,1> Xs = {0,0,0,0,0,0};
Matrix<3,6> K = {  0.435003,    0.00000,    0.0000,   0.1521001,    0.0000,    0.0000,
                  -0.000000,    0.435003,    0.0000,   0.000000,   0.1521001,    0.0000,
                   0.000000,   -0.00000,    24.80,   -0.00000,    0.0000,   0.471857}; 

                   
             //      Matrix<3,6> K = {  0.40003,    0.0000,    0.0000,   -0.28001,    0.0000,    0.0000,
             //     -0.0000,    0.2085,    0.0000,   0.0000,    -0.1061,    0.0000,
              //     0.0000,   -0.0000,    24.80,   -0.0000,    0.0000,   0.471857}; 
/////////////////////// ATTITUDE ONLY START ///////////////////////

/////////////////////// HOVER ONLY START ///////////////////////
Matrix<4,1> U_hov = {0,0,0,0}; // Output vector
Matrix<8,1> error_hov {0,0,0,0,0,0,0,0}; // State error vector
Matrix<8,1> REF_hov = {0,0,0,0,0,0,0,0}; 
Matrix<8,1> Xs_hov = {0,0,0,0,0,0,0,0};
Matrix<4,8> K_hov = {  -0.5164,   -0.0000,    0.0000,   -0.1440,   -0.0000,    0.0000,    0.0000,    0.0000,
                        0.0000,   -0.5164,    0.0000,   -0.0000,   -0.1440,   -0.0000,    0.0000,    0.0000,
                       -0.0000,    0.0000,   -0.0316,   -0.0000,    0.0000,   -0.0561,    0.0000,    0.0000,
                        0.0000,    0.0000,    0.0000,   -0.0000,    0.0000,    0.0000,   75.2536,   24.2536}; 

/////////////////////// HOVER ONLY START ///////////////////////


/////////////// Estimator matrixes Start ////////////////////////
Matrix<6,6> A = {   1,  0,  0,  DT_SEC, 0,  0,
                    0,  1,  0,  0,  DT_SEC, 0,
                    0,  0,  1,  0,  0,  DT_SEC,
                    0,  0,  0,  1,  0,  0,
                    0,  0,  0,  0,  1,  0, 
                    0,  0,  0,  0,  0,  1 };

Matrix<6,3> B = {   0.5*pow(DT_SEC,2),  0,                  0,         
                    0,                  0.5*pow(DT_SEC,2),  0,         
                    0,                  0,                  0.5*pow(DT_SEC,2),  
                    DT_SEC,             0,                  0,         
                    0,                  DT_SEC,             0,         
                    0,                  0,                  DT_SEC };
 
Matrix<6,6> H = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; 

// State vector
Matrix<6,1> Xest = {0,0,0,0,0,0};

// Prediction vector
Matrix<6,1> Xpre = {0,0,0,0,0,0};

// Measurement vector
Matrix<6,1> Z = {0,0,0,0,0,0};

// Input vector
Matrix<3,1> Uest = {0,0,0};

// Estimator gain
/* 
Matrix<6,6> Kf = {  0.0001,    0.0000,    0.0000,    0.0095,   -0.0000,    0.0000,
                    0.0000,    0.0001,   -0.0000,    0.0000,    0.0095,   -0.0000,
                    0.0000,   -0.0000,    0.1547,    0.0000,   -0.0000,    0.0000,
                    0.0000,    0.0000,    0.0000,    0.1057,   -0.0000,    0.0000,
                   -0.0000,    0.0000,   -0.0000,   -0.0000,    0.1057,   -0.0000,
                    0.0000,   -0.0000,    1.3001,    0.0000,   -0.0000,    0.0000  }; 
                    */

/* 
Matrix<6,6> Kf = {  0.0345,    0.0000,   0.0000,    0.0019,    0.0000,    0.0000,
                    0.0000,    0.0345,   0.0000,    0.0000,    0.0019,    0.0000,
                    0.0000,    0.0000,   0.0274,    0.0000,    0.0000,    0.0001,
                    0.1495,    0.0000,   0.0000,    0.0216,    0.0000,    0.0000,
                    0.0000,    0.1495,   0.0000,    0.0000,    0.0216,    0.0000,
                    0.0000,    0.0000,   0.0767,    0.0000,    0.0000,    0.0004    }; 
                    */

///*
Matrix<6,6> Kf = {  0.618520, 0.000000, 0.000000, 0.000330, 0.000000, 0.000000,
                    0.000000, 0.618520, 0.000000, 0.000000, 0.000330, 0.000000,
                    0.000000, 0.000000, 0.318880, 0.000000, 0.000000, 0.000420,
                    0.162570, 0.000000, 0.000000, 0.131210, 0.000000, 0.000000,
                    0.000000, 0.162570, 0.000000, 0.000000, 0.131210, 0.000000,
                    0.000000, 0.000000, 4.190280, 0.000000, 0.000000, 0.045440   };
               //     */
/////////////// Estimator matrixes End ////////////////////////


Servo sx;                                         // x-axis servo object (Roll)
Servo sy;                                         // y-axis servo object (Pitch)
Servo edf;                                        // EDF motor ESC object
Servo rw;                                         // Reaction Wheel ESC object
LIDARLite_v3HP myLidarLite;                       // Garmin v3HP LIDAR object
Adafruit_BNO055 bno = Adafruit_BNO055(55, 0x28);  // BNO055 9-axis IMU object 

data_prev_t pdata;                                // Previous Data struct object
data_current_t data;                              // Current Data struct object
data_estimator_t estimate;                        // Current Estimator Data Object

uint8_t sys_status, self_test, sys_err; 

float curr_time, prev_time;
float previousTime, currentTime, elapsedTime;
float control_timer{0};
float sensor_timer{0};

////////////////////
float thrust_mag;

//lidar stuff
uint8_t lidarLiteAdd{0x62};
uint16_t distance;
uint8_t lidar_measurement_response;
float temp_vec_rotated[2];

float T[2], torque[2];
float servoang1, Tmag;
float servoang2 ;
float mst{0.0f};                                  // mst: Mission Start Time
bool startFlag{false};
float tavastime{0.00f};
float estTime{0.00f};
//float p[3] = {0};



// function definitions (will need to add these to the header file but this is just for initial testing puroses to keep everything in one file)
float limit(float value, float min, float max);
float d2r(float deg);
float r2d(float rad);
void control_attitude(float r, float p, float y, float gx, float gy, float gz);
void control_attitude_hov(float r, float p, float y, float gx, float gy, float gz, float z, float vz);
float LPF( float new_sample, float old_sample, float prev_output );
float averaging(float new_sample, float old_sample);
void writeXservo(float angle);
void writeYservo(float angle);
void init_servos(void);
void writeEDF(float Ft);
float ft2omega(float Ft);
int omega2pwm(float omega);
void emergency_check(float r, float p);
void init_edf(void);
float servoRateLimit(float new_sample, float old_sample);
float IIR( float newSample, float prevOutput, float alpha);
float deadband(float new_sample, float old_sample);
void suspend(void);
void samplebno(void);
void initBno(void);
void printSerial(void);
void initLidar(void);
void sampleLidar(void);
uint8_t distanceFast(uint16_t * distance);
void rotate_to_world( float * vector );
void printLoopTime(void);
uint8_t distanceContinuous(uint16_t * distance);
void run_kalman_estimator();
void run_estimator();



void setup(void)
{
  Serial.begin(115200);

  // Initialize all hardware modules 
  initBno();
  initLidar();
  init_servos();
  //init_edf();
}


void loop(void)
{
  //printLoopTime();

  //run this statement once when loop start to grab the start time with millis()
  if(startFlag == false)
  {
    mst = millis();   
    startFlag = true;
  }

  //run BNO as fast as the loop can
  //samplebno();
  
  //SENSOR TIMER
  if(millis() - sensor_timer >= DT_MSEC)
  {
    sensor_timer = millis();
    samplebno();
    sampleLidar();
    //float temp_vec_rotated[2];
      
  }
    if(millis() - estTime >= DT_MSEC)
  {
    estTime = millis();
    run_estimator();
  }

 // CONTROLLER TIMER
  if(millis() - control_timer >= DT_MSEC)
  {
    control_timer = millis(); 
    REF_hov = {0,0,0, 0,0,0, 1.00f,0};
    //control_attitude(data.roll, data.pitch, data.yaw, data.gx, data.gy, data.gz);
    control_attitude_hov(data.roll, data.pitch, data.yaw, data.gx, data.gy, data.gz, data.z, data.vzl);
  }


  // RUN THROUGH STEP RESPONSES. 
  /*
  if(millis() - mst >= 0000 && millis() - mst <= 5000) REF = {0.00f, d2r(0.00f), 0.00f, 0.00f, 0.00f, 0.00f};
  if(millis() - mst >= 5000 && millis() - mst <= 8000) REF = {0.00f, d2r(3.00f), 0.00f, 0.00f, 0.00f, 0.00f};
  if(millis() - mst >= 8000 && millis() - mst <= 11000) REF = {d2r(-3.00f), 0.00f ,  0.00f, 0.00f, 0.00f, 0.00f};
  if(millis() - mst >= 11000 && millis() - mst <= 15000) REF = {d2r(3.00f), d2r(3.00f), 0.00f, 0.00f, 0.00f, 0.00f};

  //stop testing after 10 seconds to save battery life. 
  if(millis() - mst > 15000) suspend();
  */

 //PRINT/LIDAR TIMER
  if(millis() - tavastime >= DT_MSEC*4)
  {
    tavastime = millis();
    //sampleLidar();
    printSerial();
  }

}

void printLoopTime(void)
{
  prev_time = curr_time;
  curr_time = millis();
  float dt = (curr_time - prev_time) / 1000;
  Serial.print(dt,7); 
  Serial.println("\t");

}

void initBno(void)
{
  Serial.println("Orientation Sensor Test"); Serial.println("");
  bno.begin(OPERATION_MODE_NDOF);
  Wire.setClock(400000UL); // Set I2C frequency to 400kHz (for Arduino Due)
  bno.setExtCrystalUse(false);

  //Load the experimentally gathered calibration values into the offset struct
  adafruit_bno055_offsets_t bnoOffset;
  /*
  ////////////////BNO onboard Hopper//////////////
  //acceleration 
  bnoOffset.accel_offset_x = -20;
  bnoOffset.accel_offset_y = -176;
  bnoOffset.accel_offset_z = -33;
  //mag 
  bnoOffset.mag_offset_x = 6465;
  bnoOffset.mag_offset_y = -2144;
  bnoOffset.mag_offset_z = -724;
  //gyro
  bnoOffset.gyro_offset_x = -1;
  bnoOffset.gyro_offset_y = -1;
  bnoOffset.gyro_offset_z = -1;
  //radii
  bnoOffset.accel_radius = -1000;
  bnoOffset.mag_radius = -776;
  */
////////////////BNO onboard Hopper//////////////

///*
////////////////BNO offboard Hopper//////////////
  //acceleration 
  bnoOffset.accel_offset_x = 22;
  bnoOffset.accel_offset_y = -76;
  bnoOffset.accel_offset_z = -4;
  //mag 
  bnoOffset.mag_offset_x = 110;
  bnoOffset.mag_offset_y = 89;
  bnoOffset.mag_offset_z = -1069;
  //gyro
  bnoOffset.gyro_offset_x = -1;
  bnoOffset.gyro_offset_y = -1;
  bnoOffset.gyro_offset_z = 0;
  //radii
  bnoOffset.accel_radius = 1000;
  bnoOffset.mag_radius = 485;
////////////////BNO OFFboard Hopper//////////////
  delay(20);
//  */

  /* Initialise the sensor */
  if (!bno.begin(OPERATION_MODE_NDOF))
  {
    /* There was a problem detecting the BNO055 ... check your connections */
    Serial.print("Ooops, no BNO055 detected ... Check your wiring or I2C ADDR!");
    //while (1);
  }
  //Run the calibration loop to attain offsets 
 // calibratebno(bnoOffset);            //Comment this out if you have experimental values. 
  //Load offset type to BNO registers 
  bno.setSensorOffsets(bnoOffset); 
  //bno.setSensorOffsets(calibratebno());             //This will use the new calib function in setup() 
  delay(500);
}

void init_servos(void)
{
  //attach servo pins 
  sx.attach(XSERVO_PIN);
  sy.attach(YSERVO_PIN);
  rw.attach(RW_PIN);
  edf.attach(EDF_PIN, 1000, 2000);
  delay(200);

  //Zero Servos 
  writeXservo(0);
  writeYservo(0);
  delay(200);
}

void init_edf(void)
{
  //initialize the edf motor and run it at 1500us for 5 seconds to initialize the govenor mode to linearize the throttle curve. 
  edf.writeMicroseconds(EDF_OFF_PWM); 
  delay(2000);

  //go to 1500 and wait 5 seconds
  edf.writeMicroseconds(EDF_MIN_PWM);
  delay(5000);
}

void printSerial(void)
{
  /*
  //Serial.print("r: ");
  Serial.print(r2d(data.roll));
  Serial.print(",");  
  Serial.print(r2d(data.pitch));
  Serial.print(",");
  Serial.print(r2d(data.yaw));
  Serial.print(",  \t");

  Serial.print(data.ax);
  Serial.print(",");
  Serial.print(data.ay);
  Serial.print(",");
  Serial.print(data.az);
  Serial.print(",  \t");

  Serial.print(data.vx);
  Serial.print(",");
  Serial.print(data.vy);
  Serial.print(",");
  Serial.print(data.vzi);
  Serial.print(",");
  Serial.print(data.vzl);
  Serial.print(",  \t");

  Serial.print(data.x);
  Serial.print(",");
  Serial.print(data.y);
  Serial.print(",");
  Serial.print(data.z);
  Serial.print(",");
  Serial.print(data.zraw);
  Serial.print(",  \t");
  */

 // /*
  Serial.print(r2d(data.roll));
  Serial.print(",");  
  Serial.print(r2d(data.pitch));
  Serial.print(",");
  Serial.print(r2d(data.yaw));
  Serial.print(",  \t");
  
  Serial.print(data.z);
  Serial.print(",");
  Serial.print(100*data.vzl);
  Serial.print(",  \t");
  
  Serial.print(r2d(U_hov(0)));
  Serial.print(",");
  Serial.print(r2d(data.ang1));
  Serial.print(",");
  Serial.print(r2d(U_hov(1)));
  Serial.print(",");
  Serial.print(r2d(data.ang2));
  Serial.print(",  \t");
  
  Serial.print(U_hov(2));
  Serial.print(",");
  Serial.print(U_hov(3));
  Serial.print(",");
  Serial.print(data.Tz);
  Serial.print(",");
  Serial.print(data.Tmag);
 // */

  // Serial.print(",       ");
  // Serial.print((double) estimate.vx);
  // Serial.print(",");
  // Serial.print((double) estimate.vy);
  // Serial.print(",");
  // Serial.print((double) estimate.vz);
  Serial.println(",");
  // Serial.print(r2d(U(1)));
  // Serial.print(",");
  // Serial.print(r2d(pdata.Roll));
  // Serial.print(",");
  // Serial.print(r2d(pdata.Pitch));
  // Serial.println(",");
  
  //Serial.print("\t");
}

void calibratebno(adafruit_bno055_offsets_t offsetType)
{
  /* Display calibration status for each sensor. */
  uint8_t syst, Gyro, accel, mag = 0;
  //adafruit_bno055_offsets_t bnoO;
  bool done_flag = false;
  uint8_t count{3};         //The number of times you want to achieve 3,3,3,3 calib status

  while(done_flag == false)
  {
    bno.getCalibration(&syst, &Gyro, &accel, &mag);
    //Serial.print("CALIBRATION: Sys=");
    bno.getSensorOffsets(offsetType);

    Serial.print(syst, DEC);
    Serial.print(", ");
    Serial.print(Gyro, DEC);
    Serial.print(", ");
    Serial.print(accel, DEC);
    Serial.print(", ");
    Serial.println(mag, DEC);

    if(syst == 3 && Gyro == 3 && accel == 3 && mag == 3){
      count--;      // decrement count 
      //once reached the count dec finish task 
      if(count == 0) done_flag = true;
    }
    delay(100);
  }
  //return bnoO; 
  Serial.println(offsetType.accel_offset_x); Serial.println(offsetType.accel_offset_y); Serial.println(offsetType.accel_offset_z);
  Serial.println(" ");
  Serial.println(" ");
  Serial.println(offsetType.mag_offset_x); Serial.println(offsetType.mag_offset_y); Serial.println(offsetType.mag_offset_z);
  Serial.println(" ");
  Serial.println(" ");
  Serial.println(offsetType.gyro_offset_x); Serial.println(offsetType.gyro_offset_y); Serial.println(offsetType.gyro_offset_z);
  Serial.println(" ");
  Serial.println(" ");
  Serial.println(offsetType.accel_radius);
  Serial.println(offsetType.mag_radius);
  Serial.println("calibration is done...");
}

//This function samples the BNO055 and calculates the Euler angles and 
//gyro values, filtered, and magnetometer values for yaw/heading 
void samplebno(void){
  
  /////////EULER ANGLE////////// 
  pdata.roll = data.roll; 
  pdata.pitch = data.pitch; 
  pdata.yaw = data.yaw; 
  imu::Quaternion q = bno.getQuat();
  q.normalize();
  //convert quaternions into readable Euler angles in RADIANS 
  //given conversion from quaternions, the angles will automatically be in radians
  imu::Vector<3> euler = q.toEuler();
  //load data struct and switch z and x axes (still not sure why this is like this)
  data.roll = isnan(euler.z()) ? pdata.roll : euler.z();
  data.pitch = isnan(euler.y()) ? pdata.pitch : euler.y();
  data.yaw = isnan(euler.x()) ? pdata.yaw : euler.x();
  // data.roll = euler.z();
  // data.pitch = euler.y();
  // data.yaw = euler.x();

  //////////GYROSCOPE///////// 
  //load previous gyro data before retreiving the current data
  pdata.gx = data.gx; 
  pdata.gy = data.gy; 
  pdata.gz = data.gz;
  //read gyro values from bno
  imu::Vector<3> gyro = bno.getVector(Adafruit_BNO055::VECTOR_GYROSCOPE);
  //load temp vecotor to IIR filter the gyro values 
  data.gx = IIR(d2r(gyro.x()), pdata.gx, 0.10);
  data.gy = IIR(d2r(gyro.y()), pdata.gy, 0.10); 
  data.gz = IIR(d2r(gyro.z()), pdata.gz, 0.55);

  //////////ACCELERATION/////////
  pdata.ax = data.ax; 
  pdata.ay = data.ay; 
  pdata.az = data.az;
  imu::Vector<3> accel = bno.getVector(Adafruit_BNO055::VECTOR_LINEARACCEL);
  data.ax = IIR(accel.x(), pdata.ax, .05);
  data.ay = IIR(accel.y(), pdata.ay, .05);
  data.az = IIR(accel.z(), pdata.az, .05);


  //integrate acell vector to get velocity vector in body frame 
//  pdata.vx = data.vx;
//  pdata.vy = data.vy;
//  pdata.vz = data.vz;
//  //performe both IIR filtering and integration in 1 line.
//  data.vx = IIR((data.ax - pdata.ax)*DT_SEC, pdata.vx, .0);
//  data.vy = IIR((data.ay - pdata.ay)*DT_SEC, pdata.vy, .0);
//  //data.vz = IIR((data.az - pdata.az)*DT_SEC, pdata.vz, .20);
  
  // float atemp[2] = {0}; 
  // atemp[0] = data.ax; 
  // atemp[1] = data.ay; 
  // atemp[2] = data.az;
  // rotate_to_world( atemp ); 
  // //save previous values to prev struct
  // pdata.axw = data.axw; 
  // pdata.ayw = data.ayw;
  // pdata.azw = data.azw;
  //load current values in to current_data_t struct
  // data.axw = atemp[0];
  // data.ayw = atemp[1];
  // data.azw = atemp[2];

  //filter and differentiate in one step. 
  //data.vx = 100 * IIR(DT_SEC * (atemp[0] - pdata.axw), pdata.vx, .25);
  //data.vy = 100 * IIR(DT_SEC * (atemp[1] - pdata.ayw), pdata.vy, .25);
  //data.vz = 100 * IIR(DT_SEC * (atemp[2] - pdata.azw), pdata.vz, .25);
  
  

}

void control_attitude(float r, float p, float y, float gx, float gy, float gz)
{
  //load state vecotr 
  Xs = {r, p, 1.00f, gx, gy, 0};

  //run controller 
  error = Xs-REF; 
  U = -K * error; 

  //load desired torque vector
  //  float tx = U(2)*sin(U(0))*COM_TO_TVC;
  //  float ty = U(2)*sin(U(1))*COM_TO_TVC;
  //  float tz = U(2);

  //load new Thrust Vector from desired torque
  float Tx{-U(2)*sin(U(0))}; 
  float Ty{-U(2)*sin(U(1)*cos(U(0)))}; 
  float Tz{U(2)*cos(U(1))*cos(U(0))};           //constant for now, should be coming from position controller 

  float Tm = sqrt(pow(Tx,2) + pow(Ty,2) + pow(Tz,2)); 
  
  //different way to deduce servo angles from body forces (got this from a paper I can send you)
  //pdata.Roll = asin(-Tx/(Tmag - pow(Ty,2)));
  //pdata.Pitch = asin(Ty/Tmag);

  U(0) = asin(Tx/(Tm));
  U(1) = asin(Ty/Tm);

  // float u1temp = averaging(U(1), prev(4));
  // float u2temp = averaging(U(2), prev(5));

  //deadband 
  //U(0) = deadband(U(0), pdata.u1);
  //U(1) = deadband(U(1), pdata.u2);

  //limit the speed of servos
  //TODO: test the change made to the maxStep value 
  //      given the correction to the servo angle to TVC angle 
  //U(0) = servoRateLimit(U(0), pdata.u1);
  //U(1) = servoRateLimit(U(1), pdata.u2);

  //filter servo angles, the more filtering, the bigger the delay 
  U(0) = IIR(U(0), pdata.u1, .05);
  U(1) = IIR(U(1), pdata.u2, .05); 

  //limit servo angles to +-15º
  U(0) = limit(U(0), d2r(-15), d2r(15));
  U(1) = limit(U(1), d2r(-15), d2r(15)); 

  //write the actuation angles to the servos 
  writeXservo(r2d(U(0)));
  writeYservo(r2d(U(1)));
  //Write the thrust value (N) to the EDF motor
  writeEDF(Tm);

  //load previous data struct for filtering etc. 

  pdata.u1 = U(0);
  pdata.u2 = U(1);
  pdata.u3 = U(2);
}

void control_attitude_hov(float r, float p, float y, float gx, float gy, float gz, float z, float vz)
{
  // load prev struct with prev values 
  pdata.u1 = data.u1;
  pdata.u2 = data.u2;
  pdata.u3 = data.u3;
  pdata.u4 = data.u4;
  //load state vecotr 
  Xs_hov = {r, p, 0, gx, gy, 0, z, vz};

  //run controller 
  error_hov = Xs_hov-REF_hov; 
  U_hov = -K_hov * error_hov; 

 

  //load desired torque vector
  //  float tx = U(2)*sin(U(0))*COM_TO_TVC;
  //  float ty = U(2)*sin(U(1))*COM_TO_TVC;
  //  float tz = U(2);

  //load new Thrust Vector from desired torque
  float Tx{ -U_hov(3) * sin(U_hov(0)) }; 
  float Ty{ -U_hov(3) * sin(U_hov(1)) * cos(U_hov(0)) }; 
  float Tz{  U_hov(3) * cos(U_hov(1)) * cos(U_hov(0)) };           //constant for now, should be coming from position controller 

  float Tm = sqrt(pow(Tx,2) + pow(Ty,2) + pow(Tz,2)); 
  thrust_mag = Tm; 
  
  //different way to deduce servo angles from body forces (got this from a paper I can send you)
  //pdata.Roll = asin(-Tx/(Tmag - pow(Ty,2)));
  //pdata.Pitch = asin(Ty/Tmag);

  //U_hov(0) = asin(Tx/(Tm));
  //U_hov(1) = asin(Ty/Tm);

  //filter servo angles, the more filtering, the bigger the delay 
  U_hov(0) = IIR(U_hov(0), pdata.u1, .05);
  U_hov(1) = IIR(U_hov(1), pdata.u2, .05); 

  //limit servo angles to +-15º
  U_hov(0) = limit(U_hov(0), d2r(-15), d2r(15));
  U_hov(1) = limit(U_hov(1), d2r(-15), d2r(15)); 

  //write the actuation angles to the servos 
  writeXservo(r2d(U_hov(0)));
  writeYservo(r2d(U_hov(1)));
  writeEDF(Tm);

  
  //load previous data struct for filtering etc. 
  data.u1 = U_hov(0);
  data.u2 = U_hov(1);
  data.u3 = U_hov(2);
  data.u4 = U_hov(3);
  data.Tz = Tz; 
  data.Tmag = Tm; 
  data.ang1 = asin(Tx/Tm); 
  data.ang2 = asin(Ty/Tm);
}

float r2d(float rad)
{
  return rad * 180.00f / PI;
}

float d2r(float deg)
{
  return deg * PI / 180.00f;
}

//limit actuation of servos 
float limit(float value, float min, float max)
{
  if(value >= max ) value = max; 
  if(value <= min ) value = min; 

  return value; 
}


//Filters. Testing different ones 
float LPF( float new_sample, float old_sample, float prev_output )
{
  return ( 0.025f * new_sample + 0.025f * old_sample + 0.95f * prev_output );  
}

float IIR( float newSample, float prevOutput, float alpha)
{
  return ( (1.0f-alpha)*newSample + alpha * prevOutput);
}

float averaging(float new_sample, float old_sample)
{
  return ((new_sample + old_sample)/2.0f);
}

//Write to X servo 
void writeXservo(float angle)
{
  //map angle in degrees to pwm value for servo 
  int pwmX{ round( ( angle * X_P1 ) + X_P2 ) }; 
  sx.writeMicroseconds(pwmX); 
}

//write to Y servo 
void writeYservo(float angle)
{
  //map angle in degrees to pwm value for servo 
  int pwmY{ round( ( angle * Y_P1 ) + Y_P2 ) }; // using polynomial regression coefficients to map tvc angle to pwm vals
  sy.writeMicroseconds(pwmY); 
}

//to write command to EDF ESC
void writeEDF(float Ft)
{
  float omega{(Ft - RAD2N_P2)/RAD2N_P1}; 
  int pwm{round(omega * RAD2PWM_P1 + RAD2PWM_P2)}; 
  edf.writeMicroseconds(pwm);
}

float ft2omega(float Ft)
{
  return (Ft - RAD2N_P2)/RAD2N_P1;
}

int omega2pwm(float omega)
{
  return (int) round(omega * RAD2PWM_P1 + RAD2PWM_P2);
}

// to shut everything down if we go past max set angle 
void emergency_check(float r, float p)
{
  if(r >= MAX_VEHICLE_ANGLE_DEG || r <= -MAX_VEHICLE_ANGLE_DEG || p >= MAX_VEHICLE_ANGLE_DEG || p <= -MAX_VEHICLE_ANGLE_DEG)
  {
    writeEDF(0); 
    writeXservo(0);
    writeYservo(0);
    Serial.println("Max attitude angle reached....  "); 
    Serial.println("Vehicle is in SAFE-MODE... must restart...."); 
    while(1);
  }
}


//This is to limit the speed of actuation.. by calculating how many degrees can the servo
//actuate with respect to the loop time or DT. I am getting around 3-6degrees per 0.001 second. 
//servo: 0.1s/60º  or 0.001667s/1º 
// new_sample, old_sample are in radians. 
float servoRateLimit(float new_sample, float old_sample)
{
  // maxSteps = loop interval / (time/(time/1˚ of actuation))
  // maxSteps = dt / 0.001667  SERVO_MAX_SEC_PER_DEG = 0.001667
  // so when we do dt/0.001667 we get the degrees we can actuate 
  //dt must be in seconds 
  float maxSteps{(DT_MSEC/1000.00f)/(SERVO_MAX_SEC_PER_RAD/3.00f)};
  float temp{new_sample};

  if(new_sample >= old_sample + maxSteps) temp = old_sample + maxSteps;
  if(new_sample <= old_sample - maxSteps) temp = old_sample - maxSteps; 

  return temp; 
}

//new_sample and prev_sample are in radians 
float deadband(float new_sample, float old_sample)
{
  float db_angle{d2r(0.10f)};      //deadband angle is 0.1 degrees 
  float temp{0};
  
  if(abs(old_sample - new_sample) <= db_angle ) temp = new_sample; 
  
  return temp;

}

void suspend(void)
{
  //turn off EDF motor
  writeEDF(0); 
  //zero TVC actuators.
  writeXservo(0);
  writeYservo(0);
  //Serial.println("Max attitude angle reached....  "); 
  //Serial.println("Vehicle is in SAFE-MODE... must restart...."); 
  while(1);

}

//////////////////////////////////////////////////////////////////////////
// the functions below will be added to the sensors.cpp/sensors.h class 
// just testing initial operations and integration before getting into software
// writing. 

void initLidar(void)
{
  Wire.begin();
  // "2" = short-range, fast speed 
  myLidarLite.configure(2, lidarLiteAdd); 
}

void sampleLidar(void)
{
  lidar_measurement_response = distanceContinuous(&distance);
  if(lidar_measurement_response == 1) 
  {
    pdata.zb = data.zraw;
    data.zraw = distance/100.00f;    // uncomment for body measurement
    //data.zraw = data.z; 
//    float p[2] = {0};
//    p[2] = (float)  distance / 100.00f;
//    rotate_to_world( p );
//    data.zraw = p[2]; 
    //pdata.vz = data.vz; 
   // data.vz = IIR((data.z - pdata.z)  / DT_SEC, pdata.vz, .10);
  }

  
}

uint8_t distanceFast(uint16_t * distance)
{
  // 1. Wait for busyFlag to indicate device is idle. This must be
  //    done before triggering a range measurement.
  myLidarLite.waitForBusy();

  // 2. Trigger range measurement.
  myLidarLite.takeRange();
  myLidarLite.waitForBusy();

  // 3. Read previous distance data from device registers.
  //    After starting a measurement we can immediately read previous
  //    distance measurement while the current range acquisition is
  //    ongoing. This distance data is valid until the next
  //    measurement finishes. The I2C transaction finishes before new
  //    distance measurement data is acquired.
  *distance = myLidarLite.readDistance();

  return 1;
}

uint8_t distanceContinuous(uint16_t * distance)
{
    uint8_t newDistance = 0;

    // Check on busyFlag to indicate if device is idle
    // (meaning = it finished the previously triggered measurement)
    if (myLidarLite.getBusyFlag() == 0)
    {
        // Trigger the next range measurement
        myLidarLite.takeRange();

        // Read new distance data from device registers
        *distance = myLidarLite.readDistance();

        // Report to calling function that we have new data
        newDistance = 1;
    }

    return newDistance;
}

uint8_t distanceSingle(uint16_t * distance)
{
    // 1. Wait for busyFlag to indicate device is idle. This must be
    //    done before triggering a range measurement.
    myLidarLite.waitForBusy();

    // 2. Trigger range measurement.
    myLidarLite.takeRange();

    // 3. Wait for busyFlag to indicate device is idle. This should be
    //    done before reading the distance data that was triggered above.
    myLidarLite.waitForBusy();

    // 4. Read new distance data from device registers
    *distance = myLidarLite.readDistance();

    return 1;
}
//////////////////////////////////////////////////////////////////////////////

void rotate_to_world( float * vector )
{
  float p = data.roll;
  float q = data.pitch;
  float u = data.yaw;

  Matrix<3,1> in = { vector[0], vector[1], vector[2] };
  Matrix<3,1> out = {0,0,0};

  Matrix<3,3> R = {   cos(q)*cos(u), sin(p)*sin(q)*cos(u)-cos(p)*sin(u), cos(p)*sin(q)*cos(u)+sin(p)*sin(u) ,
                      cos(q)*sin(u), sin(p)*sin(q)*sin(u)+cos(p)*cos(u), cos(p)*sin(q)*sin(u)-sin(p)*cos(u) ,
                      -sin(q),       sin(p)*cos(q),                      cos(p)*cos(q)                      };

  out = R * in;

  vector[0] = out(0);
  vector[1] = out(1);
  vector[2] = out(2);

}

void run_kalman_estimator(){

    /* ---- Sensor processing ---- */
    float p[3] = {0}; // Position vector (z body to world)
    float v[3] = {0}; // Velocity vector (vx, vy to world)
    float a[3] = {0}; // Acceleration vector (ax, ay, az to world)

    // Rotate lidar measurement to world frame
    p[2] = data.z;
    rotate_to_world( p );

    // Perform gyrocompensation on flow and rotate to world frame.
   // v[0] = data.vx ;
  //  v[1] = data.vy ; 
  //  rotate_to_world( v );

    // Rotate acceleration to world frame
    a[0] = data.ax; a[1] = data.ay; a[2] = data.az;
    rotate_to_world( a );


    /* ---- Estimation ---- */
    // Fill input vector with acceleration
    H.Fill(0);
    Z.Fill(0);
    Uest =  {a[0], a[1], a[2]};

    // Fill measurement vector with data
    // if( data.status.pos == 1 ){
    //     H(0,0) = 1; H(1,1) = 1;
    //     Z(0) = data.x;
    //     Z(1) = data.y;
    // }

    if( lidar_measurement_response == 1){
        H(2,2) = 1;
        Z(2) = p[2]; // p[2]: z
    }

    // if( data.status.flow == 1){
      //  H(3,3) = 1; H(4,4) = 1;
      //  Z(3) = v[0]; // vx
     //   Z(4) = v[1]; // vy
    // }

    // Prediction, based on previous state and current input
    Xpre = A*Xest + B*Uest; 
    //Ppre = A*Pk_1*transpose(A) + Qf; 

    // Update prediction with update using measurements 
    //Kf = Ppre*transpose(H)*(H*Pk*transpose(H) + Rf);
    Xest = Xpre + Kf*(Z - H*Xpre); 
    
    // Fill estimate struct with values (for telemetry and stuff)
    estimate.x = Xest(0);
    estimate.y = Xest(1);
    estimate.z = Xest(2);
    estimate.vx = Xest(3);
    estimate.vy = Xest(4);
    estimate.vz = Xest(5);

}

void run_estimator()
{
  float atemp[2] = {0}; 
  float ltemp[2] = {0};

  //save previous values to prev struct
  pdata.axw = data.axw; 
  pdata.ayw = data.ayw;
  pdata.azw = data.azw;

  pdata.vx = data.vx; 
  pdata.vy = data.vy; 
  pdata.vzi = data.vzi; 
  pdata.vzl = data.vzl;

  pdata.x = data.x; 
  pdata.y = data.y; 
  pdata.z = data.z;

  atemp[0] = data.ax; 
  atemp[1] = data.ay; 
  atemp[2] = data.az;
  rotate_to_world( atemp ); 

  //load current values in to current_data_t struct
  data.axw = atemp[0];
  data.ayw = atemp[1];
  data.azw = atemp[2];

  //filter and differentiate in one step. 
  data.vx =  IIR(DT_SEC * (atemp[0] - pdata.axw), pdata.vx, .05);
  data.vy =  IIR(DT_SEC * (atemp[1] - pdata.ayw), pdata.vy, .05);
  data.vzi =  IIR(DT_SEC * (atemp[2] - pdata.azw), pdata.vzi, .05);

  data.x = IIR((data.vx - pdata.vx) * DT_SEC, pdata.x, 0.05);
  data.y = IIR((data.vy - pdata.vy) * DT_SEC, pdata.y, 0.05);
  //data.x = IIR((data.vx - pdata.vx) * DT_SEC, pdata.x, 0.05);

  // rotate lidar measurement to world and integrate for vz_idar 
  ltemp[2] = data.zraw; 
  rotate_to_world( ltemp ); 
  data.z = ltemp[2]; 

  // calculate vzworld from lidar measurement 
  data.vzl = IIR((data.z - pdata.z) / DT_MSEC, pdata.vzl, 0.05); 
}