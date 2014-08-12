#ifndef _TONECORE_H_
#define _TONECORE_H_

//transmitter sampling frequency
#define FS  22500

#define COR_L_BS 180
// Repeat number for base length
#define RPT_NUM 5
// Correlation length
#define COR_L COR_L_BS*RPT_NUM
// Head and tail length
#define TONE_HTL 23

#define COR_L_AND_HL COR_L+TONE_HTL
// Tone length
#define TONE_L (2*TONE_HTL+COR_L)
// Symbol length for data type 2
#define SYM_LEN_2 7
// Symbol length for data type 3
#define SYM_LEN_3 14

// h tone
#define HT 14
// j tone
// #define JT1A 0
// #define JT1B 2
// #define JT1C 4
// #define JT1D 6
#define JT2 8
#define JT3 10

// HPF gain
#define HPF_G 2

static int h_flag[3] = {0}; //int32(zeros(1,3));

static short int rx_reg[TONE_L] = {0}; // int16(zeros(1, tone_l));

static int rx_regI_value = 0 ;
static int *rx_regI = &rx_regI_value; //int32(0);







#endif