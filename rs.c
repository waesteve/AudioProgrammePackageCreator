#include "rs.h"

const int pp [mm+1] = { 1, 1, 0, 1} ; /* specify irreducible polynomial coeffts */

const signed char ditl_tab[28] = { 0, 7, 10, 13, 16, 23, 26,
    1, 4, 11, 14, 17, 20, 27,
    2, 5, 8, 15, 18, 21, 24,
    3, 6, 9, 12, 19, 22, 25 };

int alpha_to [nn_plusone], index_of [nn_plusone], gg [nn_sub_kk_plusone] ;
int recd [nn], data [kk], bb [nn_sub_kk] ;
int elp[nn_sub_kk_plustwo][nn_sub_kk], d[nn_sub_kk_plustwo], l[nn_sub_kk_plustwo], u_lu[nn_sub_kk_plustwo], s[nn_sub_kk_plusone] ;
int root[tt], loc[tt], z[tt_plusone], err[nn], reg[tt_plusone] ;

int decode_rs()
/* assume we have received bits grouped into mm-bit symbols in recd[i],
 i=0..(nn-1),  and recd[i] is index form (ie as powers of alpha).
 We first compute the 2*tt syndromes by substituting alpha**i into rec(X) and
 evaluating, storing the syndromes in s[i], i=1..2tt (leave s[0] zero) .
 Then we use the Berlekamp iteration to find the error location polynomial
 elp[i].   If the degree of the elp is >tt, we cannot correct all the errors
 and hence just put out the information symbols uncorrected. If the degree of
 elp is <=tt, we substitute alpha**i , i=1..n into the elp to get the roots,
 hence the inverse roots, the error location numbers. If the number of errors
 located does not equal the degree of the elp, we have more than tt errors
 and cannot correct them.  Otherwise, we then solve for the error value at
 the error location and correct the error.  The procedure is that found in
 Lin and Costello. For the cases where the number of errors is known to be too
 large to correct, the information symbols as received are output (the
 advantage of systematic encoding is that hopefully some of the information
 symbols will be okay and that if we are in luck, the errors are in the
 parity part of the transmitted codeword).  Of course, these insoluble cases
 can be returned as error flags to the calling routine if desired.   */
{
    register int i,j,u,q ;
    
    int count=0, syn_error=0;
    
    /* first form the syndromes */
    for (i=1; i<=nn_sub_kk; i++)
    { s[i] = 0 ;
        for (j=0; j<nn; j++)
            if (recd[j]!=-1)
                s[i] ^= alpha_to[(recd[j]+i*j)%nn] ;      /* recd[j] in index form */
        /* convert syndrome from polynomial form to index form  */
        if (s[i]!=0)  syn_error=1 ;        /* set flag if non-zero syndrome => error */
        s[i] = index_of[s[i]] ;
    } ;
    //printf("syn_error:%d \n" ,syn_error);
    
    //showDebugInfo();
    if (syn_error)       /* if errors, try and correct */
    {
        /* compute the error location polynomial via the Berlekamp iterative algorithm,
         following the terminology of Lin and Costello :   d[u] is the 'mu'th
         discrepancy, where u='mu'+1 and 'mu' (the Greek letter!) is the step number
         ranging from -1 to 2*tt (see L&C),  l[u] is the
         degree of the elp at that step, and u_l[u] is the difference between the
         step number and the degree of the elp.
         */
        /* initialise table entries */
        d[0] = 0 ;           /* index form */
        d[1] = s[1] ;        /* index form */
        elp[0][0] = 0 ;      /* index form */
        elp[1][0] = 1 ;      /* polynomial form */
        for (i=1; i<nn_sub_kk; i++)
        { elp[0][i] = -1 ;   /* index form */
            elp[1][i] = 0 ;   /* polynomial form */
        }
        l[0] = 0 ;
        l[1] = 0 ;
        u_lu[0] = -1 ;
        u_lu[1] = 0 ;
        u = 0 ;
        
        do
        {
            u++ ;
            if (d[u]==-1)
            { l[u+1] = l[u] ;
                for (i=0; i<=l[u]; i++)
                {  elp[u+1][i] = elp[u][i] ;
                    elp[u][i] = index_of[elp[u][i]] ;
                }
            }
            else
            /* search for words with greatest u_lu[q] for which d[q]!=0 */
            { q = u-1 ;
                while ((d[q]==-1) && (q>0)) q-- ;
                /* have found first non-zero d[q]  */
                if (q>0)
                { j=q ;
                    do
                    { j-- ;
                        if ((d[j]!=-1) && (u_lu[q]<u_lu[j]))
                            q = j ;
                    }while (j>0) ;
                } ;
                
                /* have now found q such that d[u]!=0 and u_lu[q] is maximum */
                /* store degree of new elp polynomial */
                if (l[u]>l[q]+u-q)  l[u+1] = l[u] ;
                else  l[u+1] = l[q]+u-q ;
                
                /* form new elp(x) */
                for (i=0; i<nn_sub_kk; i++)    elp[u+1][i] = 0 ;
                for (i=0; i<=l[q]; i++)
                    if (elp[q][i]!=-1)
                        elp[u+1][i+u-q] = alpha_to[(d[u]+nn-d[q]+elp[q][i])%nn] ;
                for (i=0; i<=l[u]; i++)
                { elp[u+1][i] ^= elp[u][i] ;
                    elp[u][i] = index_of[elp[u][i]] ;  /*convert old elp value to index*/
                }
            }
            u_lu[u+1] = u-l[u+1] ;
            
            /* form (u+1)th discrepancy */
            if (u<nn_sub_kk)    /* no discrepancy computed on last iteration */
            {
                if (s[u+1]!=-1)
                    d[u+1] = alpha_to[s[u+1]] ;
                else
                    d[u+1] = 0 ;
                for (i=1; i<=l[u+1]; i++)
                    if ((s[u+1-i]!=-1) && (elp[u+1][i]!=0))
                        d[u+1] ^= alpha_to[(s[u+1-i]+index_of[elp[u+1][i]])%nn] ;
                d[u+1] = index_of[d[u+1]] ;    /* put d[u+1] into index form */
            }
        } while ((u<nn_sub_kk) && (l[u+1]<=tt)) ;
        
        u++ ;
        if (l[u]<=tt)         /* can correct error */
        {
            /* put elp into index form */
            for (i=0; i<=l[u]; i++)   elp[u][i] = index_of[elp[u][i]] ;
            
            /* find roots of the error location polynomial */
            for (i=1; i<=l[u]; i++)
                reg[i] = elp[u][i] ;
            count = 0 ;
            for (i=1; i<=nn; i++)
            {  q = 1 ;
                for (j=1; j<=l[u]; j++)
                    if (reg[j]!=-1)
                    { reg[j] = (reg[j]+j)%nn ;
                        q ^= alpha_to[reg[j]] ;
                    } ;
                if (!q)        /* store root and error location number indices */
                { root[count] = i;
                    loc[count] = nn-i ;
                    count++ ;
                };
            } ;
            if (count==l[u])    /* no. roots = degree of elp hence <= tt errors */
            {
                /* form polynomial z(x) */
                for (i=1; i<=l[u]; i++)        /* Z[0] = 1 always - do not need */
                { if ((s[i]!=-1) && (elp[u][i]!=-1))
                    z[i] = alpha_to[s[i]] ^ alpha_to[elp[u][i]] ;
                else if ((s[i]!=-1) && (elp[u][i]==-1))
                    z[i] = alpha_to[s[i]] ;
                else if ((s[i]==-1) && (elp[u][i]!=-1))
                    z[i] = alpha_to[elp[u][i]] ;
                else
                    z[i] = 0 ;
                    for (j=1; j<i; j++)
                        if ((s[j]!=-1) && (elp[u][i-j]!=-1))
                            z[i] ^= alpha_to[(elp[u][i-j] + s[j])%nn] ;
                    z[i] = index_of[z[i]] ;         /* put into index form */
                } ;
                
                /* evaluate errors at locations given by error location numbers loc[i] */
                for (i=0; i<nn; i++)
                { err[i] = 0 ;
                    if (recd[i]!=-1)        /* convert recd[] to polynomial form */
                        recd[i] = alpha_to[recd[i]] ;
                    else  recd[i] = 0 ;
                }
                for (i=0; i<l[u]; i++)    /* compute numerator of error term first */
                { err[loc[i]] = 1;       /* accounts for z[0] */
                    for (j=1; j<=l[u]; j++)
                        if (z[j]!=-1)
                            err[loc[i]] ^= alpha_to[(z[j]+j*root[i])%nn] ;
                    if (err[loc[i]]!=0)
                    { err[loc[i]] = index_of[err[loc[i]]] ;
                        q = 0 ;     /* form denominator of error term */
                        for (j=0; j<l[u]; j++)
                            if (j!=i)
                                q += index_of[1^alpha_to[(loc[j]+root[i])%nn]] ;
                        q = q % nn ;
                        err[loc[i]] = alpha_to[(err[loc[i]]-q+nn)%nn] ;
                        recd[loc[i]] ^= err[loc[i]] ;  /*recd[i] must be in polynomial form */
                    }
                }
                return count ;
            }
            else {   /* no. roots != degree of elp => >tt errors and cannot solve */
                for (i=0; i<nn; i++)        /* could return error flag if desired */
                    if (recd[i]!=-1)        /* convert recd[] to polynomial form */
                        recd[i] = alpha_to[recd[i]] ;
                    else  recd[i] = 0 ;     /* just output received codeword as is */
                return -1;
            }
        }
        else{         /* elp has degree has degree >tt hence cannot solve */
            for (i=0; i<nn; i++)       /* could return error flag if desired */
                if (recd[i]!=-1)        /* convert recd[] to polynomial form */
                    recd[i] = alpha_to[recd[i]] ;
                else  recd[i] = 0 ;     /* just output received codeword as is */
            return -1;
        }}
    else {      /* no non-zero syndromes => no errors: output received codeword */
        for (i=0; i<nn; i++)
            if (recd[i]!=-1)        /* convert recd[] to polynomial form */
                recd[i] = alpha_to[recd[i]] ;
            else  recd[i] = 0 ;
        
        return 0;
    }
}


void rsInit(){
    register int i,j, mask ;
    
    mask = 1 ;
    alpha_to[mm] = 0 ;
    for (i=0; i<mm; i++)
    { alpha_to[i] = mask ;
        index_of[alpha_to[i]] = i ;
        if (pp[i]!=0)
            alpha_to[mm] ^= mask ;
        mask <<= 1 ;
    }
    index_of[alpha_to[mm]] = mm ;
    mask >>= 1 ;
    for (i=mm_plusone; i<nn; i++)
    { if (alpha_to[i-1] >= mask)
        alpha_to[i] = alpha_to[mm] ^ ((alpha_to[i-1]^mask)<<1) ;
    else alpha_to[i] = alpha_to[i-1]<<1 ;
        index_of[alpha_to[i]] = i ;
    }
    index_of[0] = -1 ;
    
    
    gg[0] = 2 ;    /* primitive element alpha = 2  for GF(2**mm)  */
    gg[1] = 1 ;    /* g(x) = (X+alpha) initially */
    for (i=2; i<=nn_sub_kk; i++)
    { gg[i] = 1 ;
        for (j=i-1; j>0; j--)
            if (gg[j] != 0)  gg[j] = gg[j-1]^ alpha_to[(index_of[gg[j]]+i)%nn] ;
            else gg[j] = gg[j-1] ;
        gg[0] = alpha_to[(index_of[gg[0]]+i)%nn] ;     /* gg[0] can never be zero */
    }
    /* convert gg[] to index form for quicker encoding */
    for (i=0; i<=nn_sub_kk; i++)  gg[i] = index_of[gg[i]] ;
  
}


int possibleControlGroups[64] ;
signed char possibleControlGroupCounts[64];
signed char difControlPackageColumns[7];
int rx_data_c3[28];

void rsDecode(int * sourceData , int type, int *resultRsData ){
    register int i = 0;
    register int j = 0;
    register int status = Failure;
    
    int loopCount = 0;
    int diffCount = 0;
    signed char mostSque = 0xFF ;
    int mostSqueIndex = -1;
    
    signed char tempCount = 0x0 ;
    signed char tempGroup[7];
    
    int resultInt = 0;
    int tempAdjustData = 0;
    
    signed char isMutiMax = 0 ;
    int tempPossibleGroupInt = 0 ;

    int * tempProgramGroup ;
    
    resultRsData[0]= SUCCESS ;
    
    switch (type) {
        case TYPE_CONTROL:
            //decodeControlPackage(needDecordData);
            
            tempAdjustData = sourceData[13];
            sourceData[13] = sourceData[9]  ;
            sourceData[9] = sourceData[12];
            sourceData[12] = sourceData[8];
            sourceData[8] = sourceData[11];
            sourceData[11] = sourceData[7];
            sourceData[7] = sourceData[10];
            sourceData[10] = tempAdjustData;
            
            for( i = 0 ; i < nn ; i ++ )
            {
 //               printf(" --%d--%d--\n ",sourceData[i] , sourceData[7+i]);
                if(sourceData[i] != sourceData[7+i])
                {
                    difControlPackageColumns[diffCount] = i ;
                    diffCount ++ ;
                }
            }
            
            loopCount = (1 <<diffCount);
            if( diffCount > 6 ){
                resultRsData[0]= -1;
                return;
            }else{
                for( i = 0 ; i < loopCount ; i++ )
                {
                    for( j =0 ; j < nn ; j++ )
                    {
                        tempGroup[j] = sourceData[j];
                    }
                    for( j = 0 ; j < diffCount ; j++ )
                    {
                        if( (i &(0x01<<j)) >0 )
                        {
                            tempGroup [ difControlPackageColumns[j] ] = sourceData[ difControlPackageColumns[j]+7 ] ;
                        }
                    }
                    for (j=0; j<nn; j++)
                    {
                        recd[j] = index_of[tempGroup[nn_subone-j]] ;
                    }
                    status = decode_rs() ;
                    if( status >=0 ){

                        tempPossibleGroupInt = 0 ;
                        for( j = 0 ; j <nn ;j++ )
                        {
                            tempPossibleGroupInt += recd[nn_subone-j]<< (j<<2);
                        }
                        possibleControlGroups[i] = tempPossibleGroupInt ;
                        possibleControlGroupCounts[i] =  0x0;
                    }else{
                        possibleControlGroupCounts[i] = 0x80;
                    }
                }
                
                for( i = 0 ; i < loopCount ; i++ )
                {
                    if ( possibleControlGroupCounts [i] == 0 )
                    {
                        tempCount = 0x0 ;
                        for( j = 0 ; j < loopCount ; j++ ){
                            if( possibleControlGroups[i] == possibleControlGroups[j]){
                                possibleControlGroupCounts[j] = 0x80;
                                tempCount++;
                            }
                        }
                        possibleControlGroupCounts[i] = tempCount ;
                        if( tempCount > mostSque )
                        {
                            mostSqueIndex = i;
                            mostSque = tempCount ;
                            isMutiMax = 0 ;
                        }else if(tempCount == mostSque ){
                            isMutiMax++;
                        }
                    }
                }
                if( mostSqueIndex < 0 || isMutiMax>0 ){
                    resultRsData[0]= -2;
                    return;
                }
                
                resultInt = possibleControlGroups[mostSqueIndex];
                resultRsData[1] = (resultInt & 0x0000000F) ;
                resultRsData[2] = (resultInt & 0x000000F0)>>4 ;
                resultRsData[3] = (resultInt & 0x00000F00)>>8 ;
            }

            
            break;
        case TYPE_RPOGRAM:
            //decodeProgrammingPackage(needDecordData);
       
            for (i = 0 ; i < 28 ; i++)
            {
                rx_data_c3[i] = sourceData[ditl_tab[i]];
            }
            for( i = 0 ;  i < 4 ; i++ ){
                tempProgramGroup = rx_data_c3 + 7*i;
                for (j=0; j<7; j++)
                {
                    recd[j] = index_of[*(tempProgramGroup+ nn_subone-j)] ;
                }
                status = decode_rs() ;
                if( status < 0){
                    resultRsData[0]=  Failure;
                    break;
                }else{
                    resultRsData[ (i*3)+1] = recd[6];
                    resultRsData[ (i*3)+2] = recd[5];
                    resultRsData[ (i*3)+3] = recd[4];
                }
            }
            
            
            break;
        default:
            break;
    }
}


void encode_rs()
/* take the string of symbols in data[i], i=0..(k-1) and encode systematically
 to produce 2*tt parity symbols in bb[0]..bb[2*tt-1]
 data[] is input and bb[] is output in polynomial form.
 Encoding is done by using a feedback shift register with appropriate
 connections specified by the elements of gg[], which was generated above.
 Codeword is   c(X) = data(X)*X**(nn-kk)+ b(X)          */
{
    register int i,j ;
    int feedback ;
    
    for (i=0; i<nn-kk; i++)   bb[i] = 0 ;
    for (i=kk-1; i>=0; i--)
    {  feedback = index_of[data[i]^bb[nn-kk-1]] ;
        if (feedback != -1)
        { for (j=nn-kk-1; j>0; j--)
            if (gg[j] != -1)
                bb[j] = bb[j-1]^alpha_to[(gg[j]+feedback)%nn] ;
            else
                bb[j] = bb[j-1] ;
            bb[0] = alpha_to[(gg[0]+feedback)%nn] ;
        }
        else
        { for (j=nn-kk-1; j>0; j--)
            bb[j] = bb[j-1] ;
            bb[0] = 0 ;
        } ;
    } ;
} ;



/***
 *param kkData:需要RS编码的值
 *param ttData：RS纠错码
 */
void rsEncode(int  *kkData , int *ttData )
{
    register int i;
    
    //doInit();
    /* for known data, stick a few numbers into a zero codeword. Data is in
     polynomial form.
     */
    for(i=0; i<kk; i++){
        data[i] = 0 ;
    }
    //init kkData;
    for( i = 0 ; i < 3 ; i++ ){
        data[ i ] =*(kkData+2-i );
    }
    
    /* encode data[] to produce parity in bb[].  Data input and parity output
     is in polynomial form
     */
    encode_rs() ;
    
    for (i=0; i<nn-kk; i++){
        *(ttData+3-i) = bb[i];
    }
    
}




