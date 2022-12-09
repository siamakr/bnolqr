//#include <Arduino.h>
#include <Wire.h>
#include <Adafruit_Sensor.h>
#include <Adafruit_BNO055.h>
#include <utility/imumaths.h>
#include <Adafruit_SPIDevice.h>
#include <Wire.h>
#include <Servo.h>
#include <Math.h>
#include <BasicLinearAlgebra.h>

//#define PI 3.141592653589793238462643383279
#define COM_TO_TVC 0.1335
#define MASS 2.458                    //Kg
#define MAX_ANGLE_SERVO 15           //Deg
#define X 0
#define Y 1
#define Z 2
//#define SERVO_MAX_SEC_PER_DEG 0.001667      //dt/.001667 |(dt=0.1) = 6º
#define SERVO_MAX_SEC_PER_DEG 0.003333f      //dt/.003333 |(dt=0.1) = 3º
#define SERVO_MAX_SEC_PER_RAD 0.0095496f      //dt/.003333 |(dt=0.1) = 3º
#define DT .01
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

using namespace BLA;

float curr_time, prev_time, dt;
float previousTime, currentTime, elapsedTime;
float control_timer{0};
float sensor_timer{0};

Matrix<3,1> U = {0,0,0}; // Output vector
Matrix<6,1> error {0,0,0,0,0,0}; // State error vector



Matrix<3,6> K = {  0.45003,    0.00000,    0.0000,   0.121001,    0.0000,    0.0000,
                  -0.000000,    0.45003,    0.0000,   0.000000,   0.121001,    0.0000,
                   0.000000,   -0.00000,    25.00,   -0.00000,    0.0000,   0.471857}; 

                   
             //      Matrix<3,6> K = {  0.40003,    0.0000,    0.0000,   -0.28001,    0.0000,    0.0000,
             //     -0.0000,    0.2085,    0.0000,   0.0000,    -0.1061,    0.0000,
              //     0.0000,   -0.0000,    24.80,   -0.0000,    0.0000,   0.471857}; 
Matrix<6,1> REF = {0,0,0,0,0,0}; 
Matrix<6,1> Xs = {0,0,0,0,0,0};

Servo sx; 
Servo sy; 
Servo edf;

typedef struct
{
float Roll, Pitch, Yaw; 
float Gx, Gy, Gz;
float Gxf, Gyf, Gzf; 
float u1, u2, u3, u4; 
} data_prev_t;

typedef struct 
{
  float Roll, Pitch, Yaw;
  float Gx, Gy, Gz;
  float ax, ay, az; 
  float vx, vy, vz; 
  float x, y, z; 
}data_current_t;


float T[2], torque[2];
  float servoang1, Tmag;
  float servoang2 ;
  int pwmx, pwmy;
  float mst{0.0f};
  bool startFlag{false};

  //define the calibration struct for BNO055



// function definitions (will need to add these to the header file but this is just for initial testing puroses to keep everything in one file)
float limit(float value, float min, float max);
float d2r(float deg);
float r2d(float rad);
void control_attitude(float r, float p, float y, float gx, float gy, float gz);
void control_attitude_red(float roll, float pitch, float yaw, float gx, float gy, float gz);
//void update_IMU(void);
//void IMU_init(void);
void calculate_IMU_error(void);
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
void printSerial(void);


/* Set the delay between fresh samples */
uint16_t BNO055_SAMPLERATE_DELAY_MS = 50;

Adafruit_BNO055 bno = Adafruit_BNO055(55, 0x28);

imu::Vector<3> euler ;
imu::Vector<3> gyro ;

data_prev_t pdata;
imu::Vector<3> ef_r;
 imu::Vector<3> gf_r;
  imu::Vector<3> gflpr_r;

  


float tavastime{0.00f};

void setup(void)
{
  Serial.begin(115200);
  Serial.println("Orientation Sensor Test"); Serial.println("");
  bno.begin(OPERATION_MODE_NDOF);
  //bno.setExtCrystalUse(true);

  //Load the experimentally gathered calibration values into the offset struct
  adafruit_bno055_offsets_t bnoOffset;
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

  delay(20);
  
  /* Initialise the sensor */
  if (!bno.begin(OPERATION_MODE_NDOF))
  {
    /* There was a problem detecting the BNO055 ... check your connections */
    Serial.print("Ooops, no BNO055 detected ... Check your wiring or I2C ADDR!");
    while (1);
  }
  //Load the offsets into the BNO055 offset registers 
  bno.setSensorOffsets(bnoOffset); 
  delay(1000);


  init_servos();
  delay(100);
  init_edf();
  //bno.setAxisRemap(Adafruit_BNO055::REMAP_CONFIG_P1);       //P1 is supposed to be default but let's see. 
  //delay(500);

}








void loop(void)
{
//      prev_time = curr_time;
//    curr_time = millis();
//    dt = (curr_time - prev_time) / 1000;
//   Serial.print(dt,7); 
//   Serial.println("\t");

  //run this statement once when loop start to grab the start time with millis()
  if(startFlag == false){
    //mst = mission start time
    mst = millis();   
    startFlag = true;
  }
  
  //interval of how fast to update servos.. in ms 
  dt = 10.00f;    //mS
 // samplebno();    // sampling BNO055 as fast as possible, might be better to sample at a constant rate
                  //see below 
  
  if(millis() - sensor_timer >= 20){
    sensor_timer = millis();
      //read BNO values and load it into imu::Vector<3> euler and imu::Vector<3> gyro 
    samplebno();
      
    
  }

  //write to servos every dt microseconds (rn set at same rate as loop so as if this if statement is not even here)
 if(millis() - control_timer >= dt){
   control_timer = millis(); 
    //write to servos 
    //load controller with new state values and run controller 
    control_attitude(ef_r.x(), ef_r.y(), ef_r.z(), gf_r.x(), gf_r.y(), gf_r.z());
    //delay(2);
  }

  // RUN THROUGH STEP RESPONSES. 
  /*
  if(millis() - mst >= 0000 && millis() - mst <= 5000) REF = {0.00f, d2r(0.00f), 0.00f, 0.00f, 0.00f, 0.00f};
  if(millis() - mst >= 5000 && millis() - mst <= 8000) REF = {0.00f, d2r(5.00f), 0.00f, 0.00f, 0.00f, 0.00f};
  if(millis() - mst >= 8000 && millis() - mst <= 11000) REF = {d2r(-5.00f), 0.00f ,  0.00f, 0.00f, 0.00f, 0.00f};
  if(millis() - mst >= 11000 && millis() - mst <= 15000) REF = {d2r(-7.00f), d2r(-7.00f), 0.00f, 0.00f, 0.00f, 0.00f};

    //stop testing after 10 seconds to save battery life. 
  if(millis() - mst > 15000) suspend();
  */
  if(millis() - tavastime >= 50){
    tavastime = millis();
    //printing roll and pitch angles of the vehicle to Serial
    //to be used in CSV file for later analysis 
    //TODO:this will need to be saved to a SD card or EEPROM when 
    //     we switch over to the TEENSY 4.2 board 
    printSerial();
      
    
  }



}

//This function samples the BNO055 and calculates the Euler angles and 
//gyro values, filtered, and magnetometer values for yaw/heading 
void samplebno(void){
  

  //EULER ANGLE 
  //load q buffer with quaternion
  imu::Quaternion q = bno.getQuat();
  q.normalize();
  //convert quaternions into readable Euler angles in RADIANS 
  //given conversion from quaternions, the angles will automatically be in radians
  imu::Vector<3> euler = q.toEuler();
  //load eulerFiltered_radians buffer with euler angles from quaternions 
  //switch x and z axes (still not sure why only the angles are switched but not the gyro)
  ef_r = {euler.z(), euler.y(), euler.x()};

  //GYROSCOPE 
  imu::Vector<3> gyro = bno.getVector(Adafruit_BNO055::VECTOR_GYROSCOPE);
  //load temp vecotor to IIR filter the gyro values 
  gf_r = {IIR(d2r(gyro.x()), pdata.Gx, 0.2), IIR(d2r(gyro.y()), pdata.Gy, 0.2), IIR(d2r(gyro.z()), pdata.Gz, 0.55)};

  //load previous gyro vectors. 
  //TODO: 
  // this needs to be changed to a rolling pointer to save memory 
  pdata.Gx = gf_r.x(); 
  pdata.Gy = gf_r.y(); 
  pdata.Gz = gf_r.z();

}

void printSerial(void){
  //Serial.print("r: ");
//  Serial.print(r2d(ef_r.x()));
//  Serial.print(",");
//  
//  Serial.println(r2d(ef_r.y()));
//  Serial.print(",");
//  Serial.print("        ");
  Serial.print(r2d(U(0)));
  Serial.print(",");
  Serial.print(r2d(U(1)));
  Serial.print(",");
  Serial.print(r2d(pdata.Roll));
  Serial.print(",");
  Serial.print(r2d(pdata.Pitch));
  Serial.println(",");
  
  //Serial.print("\t");
}

void calibratebno(void){


  /* Get the four calibration values (0..3) */
  /* Any sensor data reporting 0 should be ignored, */
  /* 3 means 'fully calibrated" */
  uint8_t system, gyro, accel, mag;
  system = gyro = accel = mag = 0;
  bno.getCalibration(&system, &gyro, &accel, &mag);

  /* The data should be ignored until the system calibration is > 0 */
  Serial.print("\t");
  if (!system)
  {
    Serial.print("! ");
  }

  /* Display the individual values */
  Serial.print("Sys:");
  Serial.print(system, DEC);
  Serial.print(" G:");
  Serial.print(gyro, DEC);
  Serial.print(" A:");
  Serial.print(accel, DEC);
  Serial.print(" M:");
  Serial.println(mag, DEC);

}
void control_attitude(float r, float p, float y, float gx, float gy, float gz){

//  Xs(0) = r;
//  Xs(1) = p;
//  Xs(2) = 1;
//  Xs(3) = gx;
//  Xs(4) = gy;
//  Xs(5) = 0;



//load state vecotr 
  Xs = {r, p, 1.00f, gx, gy, 0};

  //run controller 
  error = Xs-REF; 
  U = -K * error; 

//    Serial.print("ang1: ");
//  Serial.print(r2d(U(0)));
//  Serial.print("\t");
// Serial.print("ang2: ");
// Serial.print(r2d(U(1)));
// Serial.print("\t");

  //load desired torque vector
  //  float tx = U(2)*sin(U(0))*COM_TO_TVC;
  //  float ty = U(2)*sin(U(1))*COM_TO_TVC;
  //  float tz = U(2);

  //load new Thrust Vector from desired torque
  float Tx{-U(2)*sin(U(0))}; 
  float Ty{-U(2)*sin(U(1))}; 
  float Tz{U(2)};           //constant for now, should be coming from position controller 

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
U(0) = IIR(U(0), pdata.u1, .12);
U(1) = IIR(U(1), pdata.u2, .12); 

//limit servo angles to +-15º
U(0) = limit(U(0), d2r(-15), d2r(15));
U(1) = limit(U(1), d2r(-15), d2r(15)); 



 writeXservo(r2d(U(0)));
  //delay(.5);
 writeYservo(r2d(U(1)));
//delay(.5);
 writeEDF(Tmag);



//
// Serial.print("r:  ");
// Serial.print(r2d(r));
// Serial.print("\t");
// Serial.print("p:  ");
// Serial.print(r2d(p));
// Serial.print("\t");
// Serial.print("gx:  ");
// Serial.print(r2d(gx));
// Serial.print("\t");
// Serial.print("gy:  ");
// Serial.print(r2d(gy));
// Serial.print("\t");
// Serial.print("fgx:  ");
// Serial.print(r2d(LPF(gx, pdata.Gx, .95)));
// Serial.print("\t");
// Serial.print("fgy:  ");
// Serial.print(r2d(Xs(4)));
// Serial.print("\t");
//  Serial.print("ang1: ");
//  Serial.print(r2d(U(0)));
//  Serial.print("\t");
// Serial.print("ang2: ");
// Serial.print(r2d(U(1)));
// Serial.print("\t");
// Serial.print("Tmag: ");
// Serial.print(Tmag);
// Serial.print("\t");
// Serial.println("\t\t");

//      Serial.print("pwmx: ");
//  Serial.print(pwmx);
//  Serial.print("\t");
//     Serial.print("pwmy: ");
//  Serial.print(pwmy);
//  Serial.print("\t");
//   Serial.println("\t");

//load previous data struct for filtering etc. 
  //pdata.Roll = r; 
  //pdata.Pitch = p; 
  pdata.Yaw = y;
  pdata.Gxf = gyro.z();
  pdata.Gyf = gyro.y();
  pdata.Gzf = gyro.x();
  pdata.u1 = U(0);
  pdata.u2 = U(1);
  pdata.u3 = U(2);



}
void control_attitude_red(float roll, float pitch, float yaw, float gx, float gy, float gz){


}

float r2d(float rad){

  return rad * 180.00f / PI;
}

float d2r(float deg){

  return deg * PI / 180.00f;
}

//limit actuation of servos 
float limit(float value, float min, float max){

  if(value >= max ) value = max; 
  if(value <= min ) value = min; 

  return value; 
}
void printEvent(sensors_event_t* event) {

}


//Filters. Testing different ones 
float LPF( float new_sample, float old_sample, float prev_output ){
    return ( 0.025f * new_sample + 0.025f * old_sample + 0.95f * prev_output );  
}

float IIR( float newSample, float prevOutput, float alpha){
    return ( (1.0f-alpha)*newSample + alpha * prevOutput);
  }

float averaging(float new_sample, float old_sample){
  return ((new_sample + old_sample)/2.0f);
}



void init_servos(void){
  sx.attach(11);
  sy.attach(10);
  edf.attach(9, 1000, 2000);
  delay(200);

  //Zero Servos 
  writeXservo(0);
  writeYservo(0);


}

void init_edf(void){

  //initialize the edf motor and run it at 1500us for 5 seconds to initialize the govenor mode to linearize the throttle curve. 
  edf.writeMicroseconds(EDF_OFF_PWM); 
  delay(2000);

  //go to 1500 and wait 5 seconds
  edf.writeMicroseconds(EDF_MIN_PWM);
  delay(5000);
}

//Write to X servo 
void writeXservo(float angle){
  //map angle in degrees to pwm value for servo 
  int pwmX{ round( ( angle * X_P1 ) + X_P2 ) }; 
  sx.writeMicroseconds(pwmX); 
  pwmx = pwmX;
//   Serial.print("pwmx: ");
//  Serial.print(pwmX);
//  Serial.print("\t");
}

//write to Y servo 
void writeYservo(float angle){
  //map angle in degrees to pwm value for servo 
   int pwmY{ round( ( angle * Y_P1 ) + Y_P2 ) }; // using polynomial regression coefficients to map tvc angle to pwm vals
   pwmy = pwmY;
  sy.writeMicroseconds(pwmY); 
//   Serial.print("pwmy: ");
//  Serial.print(pwmY);
//  Serial.print("\t");
}


//to write command to EDF ESC
void writeEDF(float Ft){
  float omega{(Ft - RAD2N_P2)/RAD2N_P1}; 
  int pwm{round(omega * RAD2PWM_P1 + RAD2PWM_P2)}; 
 
  pdata.u4 = pwm;
  edf.writeMicroseconds(pwm);

}

float ft2omega(float Ft){
  return (Ft - RAD2N_P2)/RAD2N_P1;
}

int omega2pwm(float omega){
  return (int) round(omega * RAD2PWM_P1 + RAD2PWM_P2);
}


// to shut everything down if we go past max set angle 
void emergency_check(float r, float p){
  if(r >= MAX_VEHICLE_ANGLE_DEG || r <= -MAX_VEHICLE_ANGLE_DEG || p >= MAX_VEHICLE_ANGLE_DEG || p <= -MAX_VEHICLE_ANGLE_DEG){
    writeEDF(0); 
    writeXservo(0);
    writeYservo(0);
    while(1);
    Serial.println("Max attitude angle reached....  "); 
    Serial.println("Vehicle is in SAFE-MODE... must restart...."); 
  }
}


//This is to limit the speed of actuation.. by calculating how many degrees can the servo
//actuate with respect to the loop time or DT. I am getting around 3-6degrees per 0.001 second. 
//servo: 0.1s/60º  or 0.001667s/1º 
// new_sample, old_sample are in radians. 
float servoRateLimit(float new_sample, float old_sample){
  // maxSteps = loop interval / (time/(time/1˚ of actuation))
  // maxSteps = dt / 0.001667  SERVO_MAX_SEC_PER_DEG = 0.001667
  // so when we do dt/0.001667 we get the degrees we can actuate 
  //dt must be in seconds 
  float maxSteps{(dt/1000.00f)/(SERVO_MAX_SEC_PER_RAD/3.00f)};
  float temp{new_sample};

  if(new_sample >= old_sample + maxSteps) temp = old_sample + maxSteps;
  if(new_sample <= old_sample - maxSteps) temp = old_sample - maxSteps; 

  return temp; 
}

//new_sample and prev_sample are in radians 
float deadband(float new_sample, float old_sample){
  float db_angle{d2r(0.10f)};      //deadband angle is 0.1 degrees 
  float temp{0};
  
  if(abs(old_sample - new_sample) <= db_angle ) temp = new_sample; 
  

  return temp;

}

void suspend(void){
  //turn off EDF motor
  writeEDF(0); 
  //zero TVC actuators.
  writeXservo(0);
  writeYservo(0);
  
  //Serial.println("Max attitude angle reached....  "); 
  //Serial.println("Vehicle is in SAFE-MODE... must restart...."); 

  while(1);

}