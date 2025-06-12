GPSTest project

Orientation System using GPS and Compass

This project implements a modular embedded orientation system that calculates true heading by combining GPS coordinates and compass sensor data. Using UART communication, it processes input data, computes magnetic declination, and outputs corrected direction in real time.

## Features

- Real-time GPS + compass integration
- Magnetic declination correction using WMM/IGRF
- UART-based communication
- Logs heading and position data (CSV or JSON)
- Modular codebase (GPS parser, compass reader, declination module, output/logging)

## Use Cases

- Personal navigation tools
- Robotics and autonomous vehicles
- Outdoor devices (hiking, biking, UAVs)

## Technologies

- Language: C/C++
- Protocols: UART (9600–115200 baud), NMEA
- Sensors: GPS modules (e.g., Neo-6M), Compass (e.g., HMC5883L, QMC5883)

## Getting Started

1. Connect GPS and compass modules via UART
2. Flash or run the controller program on your embedded device
3. View results via serial monitor or onboard display
4. Access logs from `/logs/` directory

## Future Plans
