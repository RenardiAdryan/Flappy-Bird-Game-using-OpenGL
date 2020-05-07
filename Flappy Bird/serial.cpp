/*AUTHOR : Renardi Adryantoro P (2210161038) CEB 2016*/
/*SERIAL COMMUNICATION*/

#include "GlobalVariable.h"

struct varSerial serial;
int flagJump;

//===========================Serial Low Level==============================================================
int openSerialPort()
{
	int fd = open(SERIAL_DEV, O_RDWR | O_NOCTTY | O_NDELAY);
	if(fd == -1) return -1;
	else
	{
		struct termios port_settings;
		tcgetattr(fd, &port_settings);
		cfsetispeed(&port_settings, SERIAL_BAUD);
		cfsetospeed(&port_settings, SERIAL_BAUD);
		// set 8 bits, no parity, no stop bits
		port_settings.c_cflag &= ~PARENB;
		port_settings.c_cflag &= ~CSTOPB;
		port_settings.c_cflag &= ~CSIZE;
		port_settings.c_cflag |= CS8;
		port_settings.c_lflag |= ~ICANON;
		tcsetattr(fd, TCSANOW, &port_settings);
		fcntl(fd, F_SETFL, O_NONBLOCK);
		tcflush(fd, TCIFLUSH); /* Discards old data in the rx buffer            */
		return(fd);
	}
}

//======================================Receive Serial From Console============================================
void pollSerialPort(int serialPort)
{	
	const int bufferSize = 10;
	static unsigned char buff[bufferSize];	
	int	n;
	
	 n = read(serialPort, buff, bufferSize);
	
	if(n > 0)
	{   
		for(int i = 0; i < n; i++)
		{
			serial.count++;
			if(serial.count < 3)
			{
				serial.serialData[serial.count] = buff[i];

				if(buff[i] == '!') serial.count = 0;
			}
			else serial.count = 0;
		}

		serial.jump = int(serial.serialData[2]);
		
		flagJump = serial.jump - 48;
	}

}


void readJoystick(){
	pollSerialPort(serial.serialPort);//Pooling serial
}


//======================================Interrupt====================================
void signalHandler(int signum)
{
	cout<<"\nInterrupt Signal Received!!"<<endl;
	
	cout<<endl;
	cout<<"  ____ _____ ___  ___  			"<<endl;
	cout<<" / ___|_   _/ _ \\|  _ \\ 		"<<endl;
	cout<<" \\___ \\ | || | | | |_) |		"<<endl;
	cout<<"  ___) || || |_| |  __/ 			"<<endl;
	cout<<" |____/ |_| \\___/|_|   			"<<endl;
	cout<<endl;
	cout<<"Stop by GLUT DESTROY!!!"<<endl;
	glutDestroyWindow(0);

	exit(signum);
}

void signalTerminated(int signum)
{
	cout<<"\nTerminated Signal Received!!"<<endl;
	cout<<endl;
	cout<<"  ____ _____ ___  ___  			"<<endl;
	cout<<" / ___|_   _/ _ \\|  _ \\ 		"<<endl;
	cout<<" \\___ \\ | || | | | |_) |		"<<endl;
	cout<<"  ___) || || |_| |  __/ 			"<<endl;
	cout<<" |____/ |_| \\___/|_|   			"<<endl;
	cout<<endl;
	cout<<"Stop by GLUT DESTROY!!!"<<endl;
	cout<<"Stop by GLUT DESTROY!!"<<endl;
	glutDestroyWindow(0);
	exit(signum);
}

