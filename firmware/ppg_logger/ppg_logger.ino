// PPG logger for Seeed XIAO ESP32-C6 + HW-827 analog pulse sensor
// Samples at 1 kHz from a hardware timer, averages 5 samples -> 200 Hz output.
// Serial output: "timestamp_us,value" (value = 12-bit ADC counts, averaged).
// Wiring: sensor S -> D0 (GPIO0), + -> 3V3, - -> GND.
const int PIN = D0;
const uint32_t PERIOD_US = 1000;
const int AVG = 5;

volatile uint32_t tickTime = 0;
volatile bool ready = false;
volatile uint32_t dropped = 0;
hw_timer_t *timer;
uint32_t sum = 0; int n = 0; uint32_t lastReport = 0;

void IRAM_ATTR onTimer() {
  if (ready) dropped++;      // main loop missed a tick
  tickTime = micros();
  ready = true;
}

void setup() {
  Serial.begin(921600);
  analogReadResolution(12);
  delay(1000);
  Serial.println("timestamp_us,value");
  timer = timerBegin(1000000);
  timerAttachInterrupt(timer, &onTimer);
  timerAlarm(timer, PERIOD_US, true, 0);
}

void loop() {
  if (ready) {
    uint32_t t = tickTime; ready = false;
    sum += analogRead(PIN); n++;
    if (n == AVG) {
      Serial.print(t); Serial.print(','); Serial.println(sum / (float)AVG, 1);
      sum = 0; n = 0;
    }
  }
  if (dropped && millis() - lastReport > 1000) {
    lastReport = millis();
    Serial.print("#dropped "); Serial.println(dropped);
  }
}
