//
//  main.c
//  AudioProgrammePackageCreator
//
//  Created by jiapeng on 14-8-11.
//  Copyright (c) 2014年 jiapeng. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "ProgramPackageCreator.h"
#include "make_wav.h"

int main(int argc, const char * argv[]) {
    
    initAudioCreator();
    
    register int i = 0 ;
    int tx_data[12] ={1,2,3 ,3,2,1 , 4,5,6,6,5,4};
    createProgramAudio(tx_data);
    for( i = 0 ; i < 16*TONE_L ; i++ ){
        printf("%d     %d\r\n",i,tx_pkt[i] );
        //save 2 file
    }
    
    FILE * stream;
    stream=fopen("/Users/jiapeng/Desktop/xcodeCreatorAudio3.pcm", "w");
    fwrite(tx_pkt, sizeof(tx_pkt[0]), 16*TONE_L, stream);
    fclose(stream);
    
     write_wav("/Users/jiapeng/Desktop/test.wav", BUF_SIZE, tx_pkt, S_RATE);
    
    return 0;
}


