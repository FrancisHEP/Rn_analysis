/************************************************* 
 * Function:       SDM
 * Description:    Minimizing likelihood func, update x_temp,y_temp
 * Calls:          x_temp, y_temp in ana1to4_BiPo.cxx; LRF.h 
 * Table Accessed: 
 * Table Updated:  
 * Input:          Step Length; Truncate Number; 
 * Output:         new x_temp, y temp and the minimum value of Likelihood
 * Return:         
 * Others:          
 * *************************************************/

void SDM(double stepLength, int truncate)
{
    double x_delay, y_delay, x_delay2, y_delay2;
    double max = -1*minLHF;
    int loop = 0;
    double LKHD = max - 0.1, LKHD10 = max, LKHD01 = max, LKHD_10 = max, LKHD0_1 = max;
    while(max!=LKHD){
        if (x_temp+stepLength>324||y_temp+stepLength-237.185>-1.7321*(x_temp+stepLength-237.185)||y_temp+stepLength-237.185>-0.5773*(x_temp+stepLength-237.185)||y_temp+stepLength>324||y_temp+stepLength-237.185>0.5773*(x_temp-stepLength+237.185)||y_temp+stepLength-237.185>1.7321*(x_temp-stepLength+237.185)||x_temp-stepLength<-324||y_temp-stepLength+237.185<-1.7321*(x_temp-stepLength+237.185)||y_temp-stepLength+237.185<-0.5773*(x_temp-stepLength+237.185)||y_temp-stepLength<-324||y_temp-stepLength+237.185<0.5773*(x_temp+stepLength-237.185)||y_temp-stepLength+237.185<1.7321*(x_temp+stepLength-237.185)){break;}
        LKHD = LHF(x_temp,y_temp);
        max = LKHD;
        LKHD10 = LHF(x_temp + stepLength,y_temp);
        LKHD01 = LHF(x_temp,y_temp + stepLength);
        LKHD_10 = LHF(x_temp - stepLength,y_temp);
        LKHD0_1 = LHF(x_temp,y_temp - stepLength);
        x_delay = x_temp;
        y_delay = y_temp;
        if (LKHD10>max){max = LKHD10;x_temp = x_temp+stepLength;}
        if (LKHD01>max){max = LKHD01;y_temp = y_temp+stepLength;}
        if (LKHD_10>max){max = LKHD_10;x_temp = x_temp-stepLength;}
        if (LKHD0_1>max){max = LKHD0_1;y_temp = y_temp-stepLength;}
        if (x_temp==x_delay2&&y_temp==y_delay2) {break;}
        x_delay2 = x_delay;
        y_delay2 = y_delay;
        loop++;
        if(loop>truncate){break;}
    }
    minLHF = -max;
    
//    if(stepLength<0.1){
//        cout<<Form("Step Length is %.2f; new X and Y are: ",stepLength)<<x_temp<<"\t"<<y_temp<<", Likelihood is: "<<minLHF<<endl;
//    }
    
    return;
}

