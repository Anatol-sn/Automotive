@startuml
title Структура компонентів проєкту орієнтації в пласких координатах(GPS + Compass + UART)

package "Апаратна частина" {
  [GPS-модуль] as GPS
  [Компас] as Compass
  [UART-інтерфейс] as UART
}

package "Програмна частина" {
  [GPSProcessor] as GPSProc
  [CompassProcessor] as CompassProc
  [MagneticDeclinationModule] as Declination
  [ResultsOutput] as Output
  [StatisticsCollector] as Stats
}

GPS --> UART : координати (NMEA)
Compass --> UART : магнітне поле

UART --> GPSProc : передача GPS-даних
UART --> CompassProc : передача компас-даних

GPSProc --> Declination : координати
CompassProc --> Declination : напрямок

Declination --> Output : напрямок з поправкою
Declination --> Stats : зберігання значень

Stats --> Output : агреговані дані

@enduml


@startuml
title Системна архітектура: Орієнтування в пласких координатах (GPS + Компас)


package "Комп'ютерна система" {
  [Компонент обчислення схилення] as Declination
  [Модуль збору статистики] as Statistics
  [Модуль виводу результатів] as Output
  [Модуль обробки GPS-даних] as GPSProc
  [Модуль обробки компаса] as CompassProc
  [UART-інтерфейс] as UART
}

package "Зовнішні датчики" {
  [GPS-модуль] as GPS
  [Компас] as Compass
}

' Зв’язки між елементами
UART --> GPS
UART --> Compass

GPS --> UART
Compass --> UART

UART --> GPSProc
UART --> CompassProc

GPSProc --> Declination
CompassProc --> Declination

Declination --> Output
Declination --> Statistics

Statistics --> Output
}

@enduml

@startuml
title Орієнтування в пласких координатах (GPS + Компас)

top to bottom direction

folder "Функціональні вимоги" {
  ' Перший стовпчик
  rectangle "FR-01\n<<requirement>>" as FR01
  note bottom of FR01
    GPS через UART\n
    Пріоритет: Високий
  end note

  rectangle "FR-02\n<<requirement>>" as FR02
  note bottom of FR02
    Орієнтація з компаса\n
    Пріоритет: Високий
  end note

  ' Другий стовпчик
  rectangle "FR-03\n<<requirement>>" as FR03
  note bottom of FR03
    Магнітне схилення\n
    Пріоритет: Високий
  end note

  rectangle "FR-04\n<<requirement>>" as FR04
  note bottom of FR04
    Вивід результатів\n
    Пріоритет: Середній
  end note

  rectangle "FR-05\n<<requirement>>" as FR05
  note bottom of FR05
    Збір статистики\n
    Пріоритет: Середній
  end note

  ' Створюємо візуальне вирівнювання по горизонталі
  FR01 -[hidden]-> FR03
  FR02 -[hidden]-> FR04
}

folder "Нефункціональні вимоги" {
  ' Перший стовпчик
  rectangle "NFR-01\n<<requirement>>" as NFR01
  note bottom of NFR01
    Оновлення ≥ 1 Гц\n
    Високий
  end note

  rectangle "NFR-02\n<<requirement>>" as NFR02
  note bottom of NFR02
    UART 9600/115200\n
    Високий
  end note

  ' Другий стовпчик
  rectangle "NFR-03\n<<requirement>>" as NFR03
  note bottom of NFR03
    Кросплатформеність\n
    Середній
  end note

  rectangle "NFR-04\n<<requirement>>" as NFR04
  note bottom of NFR04
    Автовідновлення\n
    Високий
  end note

  rectangle "NFR-05\n<<requirement>>" as NFR05
  note bottom of NFR05
    Перевірка даних\n
    Високий
  end note

  NFR01 -[hidden]-> NFR03
  NFR02 -[hidden]-> NFR04
}

@enduml
