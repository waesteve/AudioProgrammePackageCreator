//
//  ProgramPackageCreator.h
//  AudioProgrammePackageCreator
//
//  Created by jiapeng on 14-8-12.
//  Copyright (c) 2014年 jiapeng. All rights reserved.
//
#include "toneCore.h"
#define sym_len 28
//% tone number
//int tone_num = sym_len/2;
#define tone_num sym_len/2

#define MAX_AUDIO 10000 //0.8*32768


#define PI 3.1415926535

extern short tx_pkt[16*TONE_L] ;
void initAudioCreator();
int createProgramAudio(int * tx_data ) ;