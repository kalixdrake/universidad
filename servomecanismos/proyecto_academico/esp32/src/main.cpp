#include <Arduino.h>
#include <driver/mcpwm.h>

// PWM por CPU (bit-banging) para observar jitter bajo carga.
static constexpr uint8_t CPU_PWM_PIN = 25;

// PWM por periferico MCPWM (hardware), normalmente mucho mas estable.
static constexpr uint8_t MCPWM_PWM_PIN = 18;

static constexpr uint32_t PWM_FREQ_HZ = 20000;  // 20 kHz
static constexpr float PWM_DUTY_PERCENT = 30.0; // 30%

static constexpr uint32_t PERIOD_US = 1000000UL / PWM_FREQ_HZ;
static constexpr uint32_t HIGH_US = (PERIOD_US * static_cast<uint32_t>(PWM_DUTY_PERCENT)) / 100UL;
static constexpr uint32_t LOW_US = PERIOD_US - HIGH_US;

// Ajustes para hacer la carga mas agresiva y aumentar jitter.
static constexpr uint32_t LOAD_INNER_ITERS = 350000;
static constexpr uint32_t LOAD_BURSTS = 4;

void busyWaitMicros(uint32_t durationUs) {
  const uint32_t start = micros();
  while ((micros() - start) < durationUs) {
    // Busy wait intentionally: this is the CPU PWM baseline.
  }
}

void cpuPwmTask(void *parameter) {
  (void)parameter;
  pinMode(CPU_PWM_PIN, OUTPUT);

  for (;;) {
    digitalWrite(CPU_PWM_PIN, HIGH);
    busyWaitMicros(HIGH_US);

    digitalWrite(CPU_PWM_PIN, LOW);
    busyWaitMicros(LOW_US);
  }
}

void cpuLoadTask(void *parameter) {
  (void)parameter;
  volatile uint32_t acc = 0;

  for (;;) {
    // Ráfagas de carga CPU para preemptar la tarea PWM y elevar jitter.
    for (uint32_t burst = 0; burst < LOAD_BURSTS; ++burst) {
      for (uint32_t i = 1; i < LOAD_INNER_ITERS; ++i) {
        acc ^= (i * 2654435761UL) + (acc >> 1);
      }
    }
    // Pausa corta para recuperar visibilidad de la PWM por CPU.
    vTaskDelay(1);
  }
}

void setupMcpwm() {
  mcpwm_gpio_init(MCPWM_UNIT_0, MCPWM0A, MCPWM_PWM_PIN);

  mcpwm_config_t config = {};
  config.frequency = PWM_FREQ_HZ;
  config.cmpr_a = PWM_DUTY_PERCENT;
  config.cmpr_b = 0.0;
  config.counter_mode = MCPWM_UP_COUNTER;
  config.duty_mode = MCPWM_DUTY_MODE_0;
  mcpwm_init(MCPWM_UNIT_0, MCPWM_TIMER_0, &config);
}

void setup() {
  Serial.begin(115200);
  delay(300);

  pinMode(CPU_PWM_PIN, OUTPUT);
  digitalWrite(CPU_PWM_PIN, LOW);

  setupMcpwm();

  Serial.println("PWM test listo.");
  Serial.println("GPIO25: PWM por CPU (espera mas jitter)");
  Serial.println("GPIO18: PWM por MCPWM (espera menos jitter)");
  Serial.printf("f=%lu Hz, duty=%.1f%%\n", PWM_FREQ_HZ, PWM_DUTY_PERCENT);
  Serial.printf("load_iters=%lu, load_bursts=%lu, load_tasks=2\n", LOAD_INNER_ITERS, LOAD_BURSTS);
  Serial.flush();

  // Ambas tareas van en el mismo core para forzar jitter en la senal por CPU.
  xTaskCreatePinnedToCore(cpuPwmTask, "cpu_pwm", 2048, nullptr, 3, nullptr, 1);
  xTaskCreatePinnedToCore(cpuLoadTask, "cpu_load_a", 2048, nullptr, 4, nullptr, 1);
  xTaskCreatePinnedToCore(cpuLoadTask, "cpu_load_b", 2048, nullptr, 4, nullptr, 1);
}

void loop() {
  delay(1000);
}