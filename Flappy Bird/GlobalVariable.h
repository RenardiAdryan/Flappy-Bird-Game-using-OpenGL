#ifndef GLOBALVARIABLE_H_
#define GLOBALVARIABLE_H_

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <termios.h>
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <csignal>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "GL/glut.h"


using namespace std;


#define SERIAL_DEV    "/dev/serial/by-id/usb-FTDI_FT232R_USB_UART_00000000-if00-port0"
#define SERIAL_BAUD    B9600
#define SERIAL_SYS    "stty -F /dev/serial/by-id/usb-FTDI_FT232R_USB_UART_00000000-if00-port0 9600"


struct varSerial
{
	int serialPort;
	int count;
	char serialData[10];
	int jump;
};
extern struct varSerial serial;



extern void pollSerialPort(int);
extern int openSerialPort();
extern void signalHandler(int);
extern void signalTerminated(int);
extern void readJoystick();
extern int flagJump;

#endif 
