#ifndef _RS_H_
#define _RS_H_

#define _RS_H_
//
//  rs.h
//  ToneCoreMain
//
//  Created by jiapeng on 14-7-21.
//  Copyright (c) 2014年 jiapeng. All rights reserved.
//


#define mm  3           /* RS code over GF(2**4) - change to suit */
#define nn  7          /* nn=2**mm -1   length of codeword */


#define tt  2          /* number of errors that can be corrected */
#define kk  3           /* kk = nn-2*tt  */

#define mm_plusone 4
#define nn_subone 6
#define nn_plusone 8
#define nn_sub_kk 4
#define nn_sub_kk_plusone 6
#define nn_sub_kk_plustwo 6

#define tt_plusone 3


#define SUCCESS 0
#define Failure -1
//param1

#define TYPE_CONTROL 2 
#define TYPE_RPOGRAM 3 

void rsInit();
void rsDecode(int * sourceData , int type, int *resultRsData );
void rsEncode(int  *kkData , int *ttData );
#endif