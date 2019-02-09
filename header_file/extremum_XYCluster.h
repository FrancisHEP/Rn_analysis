/************************************************* 
 * Function:       extreme
 * Description:    inspect if for given PMT, it's the biggest among neighbors
 * Calls:          area[110],pmtX[110],pmtY[110]
 * Table Accessed: None
 * Table Updated:  None
 * Input:          PMT index pmtIndex
 * Output:         No 
 * Return:         bool type //for PMT pmtIndex, true for it's extreme one
 * Others:          
 * *************************************************/

bool extreme(int pmtIndex, double qS2T)
{
    double max = area[pmtIndex];
    double thresh = qS2T / 55.;
    if (thresh<10)
    {
        thresh = 10;
    }
    if(max<thresh)
    {
        return false;
    }
    int index = pmtIndex;
    for (int i = 55; i < 110; i++)
    {
        if(i == pmtIndex)
        {
            continue;
        }
        if(pmtX[pmtIndex]*pmtX[pmtIndex]+pmtY[pmtIndex]*pmtY[pmtIndex]<82000)
        {
            if((pmtX[pmtIndex]-pmtX[i])*(pmtX[pmtIndex]-pmtX[i])+(pmtY[pmtIndex]-pmtY[i])*(pmtY[pmtIndex]-pmtY[i])<7000)
            {
                if(area[pmtIndex] < area[i])
                {
                    index = i;
                }
            }
        }
        else
        {
            if(((pmtX[pmtIndex]-pmtX[i])*(pmtX[pmtIndex]-pmtX[i])+(pmtY[pmtIndex]-pmtY[i])*(pmtY[pmtIndex]-pmtY[i])<7000)||((pmtX[pmtIndex]-pmtX[i])*(pmtX[pmtIndex]-pmtX[i])+(pmtY[pmtIndex]-pmtY[i])*(pmtY[pmtIndex]-pmtY[i])<20000&&pmtX[i]*pmtX[i]+pmtY[i]*pmtY[i]>82000))
            {
                if(area[pmtIndex] < area[i])
                {
                    index = i;
                }
            }
        }
    }
    if(index != pmtIndex)
    {
        return false;
    }
    else 
    {
        return true;
    }
}
