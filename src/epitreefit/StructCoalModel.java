package epitreefit;

import cern.colt.matrix.*;
import cern.colt.list.*;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Superclass for the structured coalescent model
 * @author David
 */
public class StructCoalModel {
    
    DoubleMatrix2D F;
    DoubleMatrix2D G;
    DoubleMatrix1D Y;
    int states;
    int transCount;
    int nonTransCount;
    ArrayList<ArrayList<Integer>> transCountArray = new ArrayList<ArrayList<Integer>>();
    double currTime;
    DoubleArrayList recordedTimes = new DoubleArrayList();

    
    public void make() 
    {   
        
        /**
         * 
         * NEED TO UPDATE THIS FOR EACH MODEL!!
         * 
         */
        
    }
    

    public void updateF(DoubleMatrix1D xCurr, DoubleArrayList theta, double currTime) 
    {
        /**
         * 
         * NEED TO UPDATE THIS FOR EACH MODEL!!
         * 
         */
        
    }
    
    public void updateG(DoubleMatrix1D xCurr, DoubleArrayList theta)
    {
         /**
         * 
         * NEED TO UPDATE THIS FOR EACH MODEL!!
         * 
         */
    }
    
    public void updateY(DoubleMatrix1D xCurr) 
    {
         /**
         * 
         * NEED TO UPDATE THIS FOR EACH MODEL!!
         * 
         */
    }
    
    public double computeCoalRate(int k, int l) 
    {
        double rate = (F.getQuick(k,l) + F.getQuick(l,k)) / (Y.getQuick(k)*Y.getQuick(l));
//        double rate = 0.0;
//        if (k == 0 & l == 0) {
//            rate = F.getQuick(k,l);
//        } else {
//            rate = (F.getQuick(k,l) + F.getQuick(l,k)) / (Y.getQuick(k)*Y.getQuick(l));
//        }
        return rate; 
    }
    
    public double computeCoalRateOneDirection(int k, int l) 
    {
        double rate = F.getQuick(k,l) / (Y.getQuick(k)*Y.getQuick(l));
        return rate;
        
    }
    
    public double computeMoveRate(int k, int l, int[] lineagesByState)
    {
        
        double rate = (G.getQuick(l,k)/Y.getQuick(k)) + (F.getQuick(l,k)/Y.getQuick(k))*((Y.getQuick(l)-lineagesByState[l])/Y.getQuick(l));
        return rate;
    }
    
}
