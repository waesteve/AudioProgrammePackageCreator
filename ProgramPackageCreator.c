//
//  main.c
//  AudioProgrammePackageCreator
//
//  Created by jiapeng on 14-8-11.
//  Copyright (c) 2014年 jiapeng. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "rs.h"
//#include "toneCore.h"
#include "ProgramPackageCreator.h"

int fm[32];
float tone_env[ TONE_L ] ;
short tx_pkt[16*TONE_L] ;

int f1_head = 0 ;
int f2_head = 0 ;
//% Interleave table
const int itlv_tab[sym_len] = { 0,7,14, 21, 8, 15 ,22 ,1, 16, 23, 2 ,9, 24 ,3, 10, 17, 4, 11, 18 ,25, 12, 19, 26 ,5, 20, 27, 6, 13 };


void initAudioCreator(){
    register int i = 0;
    
    int fmIndex[32] ={ 16 ,17,18,19,20,21,22,23,24,25,26,27,28,29,30 ,31, 33,35,37,39,41,43,45,47,49, 53, 55, 57, 61, 64 ,65, 67};
    
    //frequency of 32 tones
    
    for( i = 0 ; i < 32; i++ ){
        fm[i] = fmIndex[i]*FS/COR_L_BS ;
    }
    
    float tempValue = 0 ;
    for( i = 0 ; i < TONE_HTL ; i++ ){
        tempValue = i/(float)TONE_HTL ;
        tone_env[i] = tempValue;
        tone_env[TONE_L-1-i] = tempValue;
    }
    //memset( tone_env+TONE_L , 1 ,  COR_L );
    for( i = TONE_HTL ; i <COR_L_AND_HL ; i++ ){
        tone_env[i] = 1.0f;
    }
    
    //tone_env = [(0:(tone_hl-1))/tone_hl ones(1,tone_bl) ((tone_tl-1):-1:0)/tone_tl];
    
    //Frame head frequency: first tone
    f1_head = fm[14]; // h tone
    f2_head = fm[10] ;
    
    //init f1
    for( i  = 0 ; i < TONE_L ; i++ ){
        tx_pkt[i] = (MAX_AUDIO*tone_env[i])*sin(2*PI*f1_head/(double)FS*i) ;
        //printf("%d     %d    %f \r\n",i, tx_pkt[i] ,sin(2*PI*f1_head/(double)FS*i) );
    }
    //init f2
    for( i  = 0 ; i < TONE_L ; i++ ){
        tx_pkt[i+TONE_L] = (MAX_AUDIO*tone_env[i])*sin(2*PI*f2_head/(double)FS*i) ;
    }
    rsInit();
}




int createProgramAudio(int * tx_data ) {
    

    register int i = 0 ;
    register int j = 0 ;
    
    int tx_data_c[28] ={0};
    
    //% Frame head signal
    
    int tempRsEncodeData[4];
    for(i = 0 ; i < 4 ; i++ ){
        rsEncode( tx_data +i*3 , tempRsEncodeData );
        tx_data_c[i*nn+ 0] = tx_data[i*3+0];
        tx_data_c[i*nn+ 1] = tx_data[i*3+1];
        tx_data_c[i*nn+ 2] = tx_data[i*3+2];
        tx_data_c[i*nn+ 3] = tempRsEncodeData[0];
        tx_data_c[i*nn+ 4] = tempRsEncodeData[1];
        tx_data_c[i*nn+ 5] = tempRsEncodeData[2];
        tx_data_c[i*nn+ 6] = tempRsEncodeData[3];
    }

    // % Interleaving
    int tx_data_i[sym_len];
    for( i = 0 ; i < sym_len; i++){
        tx_data_i[i] = tx_data_c[itlv_tab[i]];
    }
    
    // Modulation
    int tx_f[2][tone_num] ={0} ;
    for( i = 0 ; i < nn ; i++ ){
        tx_f[0][i<<1] =tx_data_i[0+(i<<2)]+0;
        tx_f[1][i<<1] =tx_data_i[1+(i<<2)]+16;
        tx_f[0][(i<<1)+1] =tx_data_i[2+(i<<2)]+8;
        tx_f[1][(i<<1)+1] =tx_data_i[3+(i<<2)]+24;
    }
    
    int hjHeadPartEndIndex = 2*TONE_L;
    for( i = 0 ; i < tone_num ; i++ ){
        for( j = 0 ; j < TONE_L ; j ++ ){
            tx_pkt[i* TONE_L +j +hjHeadPartEndIndex ] = MAX_AUDIO*tone_env[j]/sqrt(2)
            *(sin(2*PI*fm[tx_f[0][i]]/FS*j)+sin(2*PI*fm[tx_f[1][i]]/FS*j));
            
        }
    }

    return 0;
}
