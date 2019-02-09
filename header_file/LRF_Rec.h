/************************************************* 
 * Function:       LRF_Rec
 * Description:    update X, Y from input (X0, Y0)
 * Calls:          SDM, LHF, x_temp, y_temp, minLHF
 * Table Accessed: None
 * Table Updated:  None
 * Input:          Start Position X0 Y0
 * Output:         None 
 * Return:         None
 * Others:          
 * *************************************************/

void LRF_Rec(double x0, double y0, double stepLength)
{
x_temp = x0;
y_temp = y0;
minLHF = 0;
SDM(stepLength,672/stepLength);
SDM(stepLength/10,20);
SDM(stepLength/100,20);
SDM(stepLength/1000,20);
SDM(stepLength/10000,20);
}

