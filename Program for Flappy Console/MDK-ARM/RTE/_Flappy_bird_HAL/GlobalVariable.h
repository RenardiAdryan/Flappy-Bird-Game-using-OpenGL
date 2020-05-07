#ifndef __GLOBAL_VARIABLE_
#define __GLOBAL_VARIABLE_

#define BUF_	50
#define BUF_SIZE 70

#define digitalRead(x,y)													HAL_GPIO_ReadPin(x,y)	
#define HIGH																			1
#define LOW																				0
#define digitalWrite(x,y,val)									((val) ? (HAL_GPIO_WritePin(x,y,GPIO_PIN_SET))  : (HAL_GPIO_WritePin(x,y,GPIO_PIN_RESET)) )


extern UART_HandleTypeDef huart1;
extern UART_HandleTypeDef huart3;
extern DMA_HandleTypeDef hdma_usart1_rx;
extern DMA_HandleTypeDef hdma_usart3_rx;

extern void InitSerialIMU_DMA(void);
extern void Parsing_IMU(void);
extern void resetIMU(void);
extern int dataYaw,dataPitch,dataRoll;
extern int sampleData;
#endif