package epitreefit;

import cern.colt.list.*;

/**
 *
 * @author David
 */
public class SampleWriter {
    
    String sampleString;

    public void samplesToString(DoubleArrayList theta) 
    {
        
        int cols = theta.size();
        sampleString = "";
        for (int param = 0; param < (cols - 1); ++param) {
            sampleString = sampleString + Double.toString(theta.getQuick(param)) + ", ";
        }
        sampleString = sampleString + Double.toString(theta.getQuick(cols - 1));
    }
}
