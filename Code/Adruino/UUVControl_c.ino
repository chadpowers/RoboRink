#include <Servo.h>
#include <Stepper.h>
#include <string.h>
#include <Wire.h>
#include <Adafruit_Sensor.h>
#include <Adafruit_LSM303_U.h>
#include <Adafruit_L3GD20_U.h>
#include <Adafruit_9DOF.h>

/* Assign a unique ID to the sensors */
Adafruit_9DOF                 dof   = Adafruit_9DOF();
Adafruit_LSM303_Accel_Unified accel = Adafruit_LSM303_Accel_Unified(30301);
Adafruit_LSM303_Mag_Unified   mag   = Adafruit_LSM303_Mag_Unified(30302);

/* Update this with the correct SLP for accurate altitude measurements */
float seaLevelPressure = SENSORS_PRESSURE_SEALEVELHPA;

#define STEPS 200 // number of steps on our stepper motor

Servo esc1;
Servo esc2;
Stepper stpr = Stepper(STEPS, 4, 5, 6, 7); // H-bridge feeds pins 4-7
int inputs[3];

void initSensors()
{
  if(!accel.begin())
  {
    // There was a problem detecting the LSM303 ... check your connections 
    Serial.println(F("Ooops, no LSM303 detected ... Check your wiring!"));
    while(1);
  }
  if(!mag.begin())
  {
    // There was a problem detecting the LSM303 ... check your connections 
    Serial.println("Ooops, no LSM303 detected ... Check your wiring!");
    while(1);
  }
}

void setup() {
  // attach ESCs as servos, then open serial connection
  esc1.attach(2); // plug ESC signal receiver into PWM 2
  esc2.attach(3); // plug ESC signal receiver into PWM 3
  stpr.setSpeed(20); // rate of rotation for stepper
  initSensors();
  establishContact();
}

void loop() {
  start:
  while (Serial.available() > 0) // Don't read unless there is data
  {
    for (int i = 0; i < 3; i++) {
      String servo = Serial.readStringUntil(':');
      if (servo == "-") {
        Serial.readString();
        Serial.flush();
        Serial.print("CLOSING CONNECTION");
        Serial.end();
        delay(100);
        establishContact();
        goto start;
      } else if (servo != "") {
        //here you could check the servo number
        String pos = Serial.readStringUntil('&');
        inputs[i] = pos.toInt(); // assuming integer inputs
      }
    }
    // map inputs to useful numbers for motors and stepper (with checks)
    inputs[0] = map(constrain(inputs[0], -100, 100), -100, 100, 0, 180);
    inputs[1] = map(constrain(inputs[1], -100, 100), -100, 100, 0, 180);
    inputs[2] = constrain(inputs[2], -100, 100);
    esc1.write(inputs[0]);
    esc2.write(inputs[1]);
    stpr.step(inputs[2]);

    sensors_event_t accel_event;
    sensors_event_t mag_event;
    sensors_vec_t   orientation;

    // Calculate pitch and roll from the raw accelerometer data
    accel.getEvent(&accel_event);
    
    // send acceleration data (x, y, z)
    Serial.print("x:");
    Serial.print(accel_event.acceleration.x);
    Serial.print(",y:");
    Serial.print(accel_event.acceleration.y);
    Serial.print(",z:");
    Serial.print(accel_event.acceleration.z);
    Serial.print(",");
    
    if (dof.accelGetOrientation(&accel_event, &orientation))
    {
      // 'orientation' should have valid .roll and .pitch fields
      Serial.print(F("R:"));
      Serial.print(orientation.roll);
      Serial.print(F(",P:"));
      Serial.print(orientation.pitch);
      Serial.print(F(","));
    }
  
    // Calculate the heading using the magnetometer
    mag.getEvent(&mag_event);
    if (dof.magGetOrientation(SENSOR_AXIS_Z, &mag_event, &orientation))
    {
      // 'orientation' should have valid .heading data now
      Serial.print(F("H:"));
      Serial.print(orientation.heading);
      Serial.print(F(";\t"));
    }
  }
  delay(7);
}

void establishContact() {
  Serial.begin(9600);
  Serial.setTimeout(100);
  while (!Serial) {
    ; // wait for serial port to connect. Needed for native USB port only
  }
  
  while (Serial.available() <= 0) {
    Serial.print('A');   // send a capital A
    delay(300);
  }
}
