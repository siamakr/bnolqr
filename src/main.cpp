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
#define DT .01
#define SERVO_MAX_DEGREE_PER_DT 12
#define MAX_VEHICLE_ANGLE_DEG 35.00f
#define DEADBAND_ANGLE_DEG 0.001f

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
Matrix<3,6> K = {  1.20803,    0.0000,    0.0000,   -0.10601,    0.0000,    0.0000,
                  -0.0000,    1.2085,    0.0000,   0.0000,    -0.1061,    0.0000,
                   0.0000,   -0.0000,    24.80,   -0.0000,    0.0000,   0.471857}; 
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


/* Set the delay between fresh samples */
uint16_t BNO055_SAMPLERATE_DELAY_MS = 50;

Adafruit_BNO055 bno = Adafruit_BNO055(55, 0x28);

imu::Vector<3> euler ;
imu::Vector<3> gyro ;

data_prev_t pdata;
imu::Vector<3> ef_r;
 imu::Vector<3> gf_r;
  imu::Vector<3> gflpr_r;

  




void setup(void)
{
  Serial.begin(9600);
  Serial.println("Orientation Sensor Test"); Serial.println("");
 // bno.begin();

  /* Initialise the sensor */
  if (!bno.begin())
  {
    /* There was a problem detecting the BNO055 ... check your connections */
    Serial.print("Ooops, no BNO055 detected ... Check your wiring or I2C ADDR!");
    while (1);
  }
  delay(1000);

init_servos();
delay(100);
//init_edf();
bno.setExtCrystalUse(true);
bno.setAxisRemap(Adafruit_BNO055::REMAP_CONFIG_P1);       //P1 is supposed to be default but let's see. 
delay(500);

}








void loop(void)
{
//      prev_time = curr_time;
//    curr_time = millis();
//    dt = (curr_time - prev_time) / 1000;
//   Serial.print(dt,7); 
//   Serial.println("\t");
  
  //interval of how fast to update servos.. in microseconds 
  dt = 10.00f;    //mS

  
  if(millis() - sensor_timer >= 20){
    sensor_timer = millis();
      //read BNO values and load it into imu::Vector<3> euler and imu::Vector<3> gyro 
    euler = bno.getVector(Adafruit_BNO055::VECTOR_EULER);
    gyro = bno.getVector(Adafruit_BNO055::VECTOR_GYROSCOPE);
   // gyro.z() = -gyro.z();


    //load temp vector to convert degrees to radians and switch z and x axis (orientation of imu) 
    ef_r = {d2r(euler.z()-18.75), d2r(euler.y()), d2r(euler.x())};

    //load temp vecot to IIR filter the gyro values 
     gf_r = {IIR(d2r(gyro.x()), pdata.Gx, 0.75), IIR(d2r(gyro.y()), pdata.Gy, 0.75), IIR(d2r(gyro.x()), pdata.Gz, 0.55)};
    // load temp vector to LPF filter of gyro values (not working rn, IIR filter doing a great job so far)
   // gflpr_r = {LPF(d2r(gyro.x()), pdata.Gxf, pdata.Gx), LPF(d2r(gyro.y()), pdata.Gyf, pdata.Gy), LPF(d2r(gyro.z()), pdata.Gzf, pdata.Gz)};

//     Serial.print("gxr:  ");
// Serial.print(r2d(gf_r.x()));
// Serial.print("\t");
//      Serial.print("r:  ");
// Serial.print(r2d(ef_r.x()));
// Serial.print("\t");
//// Serial.print("gyr:  ");
//// Serial.print(r2d(gf_r.y()));
//// Serial.print("\t");Ò
//      Serial.print("gx:  ");
// Serial.print(gyro.x());
// Serial.print("\t");
//// Serial.print("gyr:  ");
//// Serial.print(r2d(gyro.y()));
//// Serial.print("\t");
//   Serial.println("\t");

      pdata.Gx = gf_r.x(); 
      pdata.Gy = gf_r.y(); 
      pdata.Gz = gf_r.z();
    
  }

  //write to servos every dt microseconds (rn set at same rate as loop so as if this if statement is not even here)
 if(millis() - control_timer >= dt){
   control_timer = millis(); 
    //write to servos 
    //load controller with new state values and run controller 
    control_attitude(ef_r.x(), ef_r.y(), ef_r.z(), gf_r.x(), gf_r.y(), gf_r.z());
    //delay(2);
  }

//writeEDF(pdata.u3);
//emergency_check(euler.z(), euler.y());
//delay(BNO055_SAMPLERATE_DELAY_MS);
}
void updatebno(void){
//will put al BNO stuff here later 
}

void calibratebno(void){


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

//   Serial.print("r:  ");
// Serial.print(r2d(Xs(0)));
// Serial.print("\t");
//  Serial.print("p:  ");
// Serial.print(r2d(Xs(1)));
// Serial.print("\t");
//  Serial.print("gx:  ");
// Serial.print(r2d(Xs(3)));
// Serial.print("\t");
//  Serial.print("r:  ");
// Serial.print(r2d(Xs(4)));
// Serial.print("\t");

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
  float Tx{U(2)*sin(U(0))}; 
  float Ty{-U(2)*sin(U(1))*cos(U(0))}; 
  float Tz{U(2)};           //constant for now, should be coming from position controller 

 Tmag = sqrt(pow(Tx,2) + pow(Ty,2) + pow(Tz,2)); 
  
  //different way to deduce servo angles from body forces (got this from a paper I can send you)
   //U(0) = asin(-Tx/(Tmag - pow(Ty,2)));
   //U(1) = asin(Ty/Tmag);

// float u1temp = averaging(U(1), prev(4));
// float u2temp = averaging(U(2), prev(5));

//deadband 
//U(0) = deadband(U(0), pdata.u1);
//U(1) = deadband(U(1), pdata.u2);

//limit the speed of servos
U(0) = servoRateLimit(U(0), pdata.u1);
U(1) = servoRateLimit(U(1), pdata.u2);

//filter servo angles, the more filtering, the bigger the delay 
U(0) = IIR(U(0), pdata.u1, .20);
U(1) = IIR(U(1), pdata.u2, .20); 

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
  pdata.Roll = r; 
  pdata.Pitch = p; 
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
float servoRateLimit(float new_sample, float old_sample){
  float maxSteps{(dt)/SERVO_MAX_SEC_PER_DEG};
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