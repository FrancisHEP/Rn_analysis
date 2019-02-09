/************************************************* 
 * Function:       LHF
 * Description:    Calculate likelihood function for given x-y and LRF of 55PMTs
 * Calls:          Globle Variables in ana1to4_BiPo.cxx 
 * Table Accessed: None
 * Table Updated:  None
 * Input:          Current X-Y Position 
 * Output:         No 
 * Return:         lhf //Likelihood value for the input x0 y0. 
 * Others: 
 * Modification: 2017-04-23 disable 11201  (ii = 80)       
 * Modification: 2017-06-18 try with disable 10900 (ii = 55)        
 * *************************************************/

double LHF(double x0, double y0)
{
    double lhf = 0;
    double r_i = 0;
    double r_ii = 0;
    double r_middle = 0;
    //double eta[110] = {0};
    double P = 0;
    double ratio = 0;
    double peak = 0;
    double qS2T = 0;
    double rho0 = 0;
    double rho_image = 1;
    double rho[110] = {0};
    rho_image = 2;	//reduce ratio of the weight on A[i] of the image PMT
    rho0 = 1.525;
    P = 0;
    for (int i = 55; i < 110; i++)
    {
    //reduce ratio of the weight on A[i] of PMTs depend on sqare of their r2
    rho[i] = rho0*pow((pmtX[i]*pmtX[i]+pmtY[i]*pmtY[i])/86349.3,2);
    }
    for(int i = 55; i < 110; i++){
        //if ( i == 55 || i == 80 )
        if ( i == 80 )
        {
            continue;
        }
    if(pmtX[i]*pmtX[i]+pmtY[i]*pmtY[i]>59000){
        if(pmtX[i]*pmtX[i]+pmtY[i]*pmtY[i]>86000){ratio = (sqrt(pmtX[i]*pmtX[i]+pmtY[i]*pmtY[i])+2*(sqrt(324*324+86.8157*86.8157)-sqrt(pmtX[i]*pmtX[i]+pmtY[i]*pmtY[i])))/sqrt(pmtX[i]*pmtX[i]+pmtY[i]*pmtY[i]);}
        else{ratio = (sqrt(pmtX[i]*pmtX[i]+pmtY[i]*pmtY[i])+2*(324-sqrt(pmtX[i]*pmtX[i]+pmtY[i]*pmtY[i])))/sqrt(pmtX[i]*pmtX[i]+pmtY[i]*pmtY[i]);}
        r_i = sqrt((x0-pmtX[i])*(x0-pmtX[i])+(y0-pmtY[i])*(y0-pmtY[i]));
        r_ii = sqrt((x0-pmtX[i]*ratio)*(x0-pmtX[i]*ratio)+(y0-pmtY[i]*ratio)*(y0-pmtY[i]*ratio));
        r_middle = (ratio-1)/2*sqrt(pmtX[i]*pmtX[i]+pmtY[i]*pmtY[i]);
        peak = 2*A[i]*exp(-a[i]*(r_middle/r[i])/(1+pow((r_middle/r[i]),(1-Alpha[i])))-b[i]/(1+pow((r_middle/r[i]),(-Alpha[i]))));
        if(peak>A[i]&&A[i]*A[i]/peak<0.25){
        eta[i] = A[i]*(1-rho[i]*A[i])*exp(-a[i]*(r_i/r[i])/(1+pow((r_i/r[i]),(1-Alpha[i])))-b[i]/(1+pow((r_i/r[i]),(-Alpha[i]))))+A[i]*(1-rho_image*A[i])*exp(-a[i]*(r_ii/r[i])/(1+pow((r_ii/r[i]),(1-Alpha[i])))-b[i]/(1+pow((r_ii/r[i]),(-Alpha[i]))));
        }
        else if(peak>A[i]&&A[i]*A[i]/peak>=0.25){
        eta[i] = A[i]*(1-rho[i]*A[i])*exp(-a[i]*(r_i/r[i])/(1+pow((r_i/r[i]),(1-Alpha[i])))-b[i]/(1+pow((r_i/r[i]),(-Alpha[i]))))+A[i]*(1-rho_image*A[i])*exp(-a[i]*(r_ii/r[i])/(1+pow((r_ii/r[i]),(1-Alpha[i])))-b[i]/(1+pow((r_ii/r[i]),(-Alpha[i]))));
        }
        else{eta[i] = A[i]*(1-rho[i]*A[i])*exp(-a[i]*(r_i/r[i])/(1+pow((r_i/r[i]),(1-Alpha[i])))-b[i]/(1+pow((r_i/r[i]),(-Alpha[i]))))+A[i]*(1-rho_image*A[i])*exp(-a[i]*(r_ii/r[i])/(1+pow((r_ii/r[i]),(1-Alpha[i])))-b[i]/(1+pow((r_ii/r[i]),(-Alpha[i]))));
        }
    }
    else{
        r_i = sqrt((x0-pmtX[i])*(x0-pmtX[i])+(y0-pmtY[i])*(y0-pmtY[i]));
        r_ii = 0;
        r_middle = 0;
        peak = 1;
        eta[i] = A[i]*(1-rho[i]*A[i])*exp(-a[i]*(r_i/r[i])/(1+pow((r_i/r[i]),(1-Alpha[i])))-b[i]/(1+pow((r_i/r[i]),(-Alpha[i]))));
        }
        P += eta[i];
        qS2T += area[i];
    }
    for(int j = 55; j < 110; j++){
        //if ( j == 55 || j == 80 )
        if ( j == 80 )
        {
           continue;
        }
        lhf += area[j]/qS2T*log(eta[j]/P);
//        cout << "PMTid = " << j << '\t' << "lhf = " << area[j]/qS2T*log(eta[j]/P) << endl;
    }
    return lhf;
}

