package epitreefit;

import cern.colt.list.*;
//import java.util.ArrayList;
import cern.jet.random.*;
import cern.colt.matrix.*;

/**
 * Specifies the epidemiological model
 * Using this model for directly transmitted dengue model in the MBE paper
 * @author David
 */
public class EpiModel2PopDengue extends EpiModel {

    //DoubleArrayList params = new DoubleArrayList();
    //ArrayList<Integer> estParams = new ArrayList<Integer>(); //indexes of the estimated params
    //ArrayList<Integer> fixedParams = new ArrayList<Integer>(); //indexes of the fixed params
    DoubleArrayList currInits;
    Poisson poissDist;
    //DoubleFactory1D factory1D = DoubleFactory1D.dense;
    //DoubleFactory2D factory2D = DoubleFactory2D.dense;
    
    AbstractContinousDistribution normDist; cern.jet.random.engine.RandomEngine engine;
    
    double seasNowFocal; double seasNowGlobal; double seasHeightFocal; double seasHeightGlobal; double betaSeasFF; double betaSeasGF; double betaSeasGG; double betaSeasFG;
    double betaFF; double betaGF; double betaGG; double betaFG; double alphaFocal; double alphaGlobal; double deltaFocal; double deltaGlobal;
    double NFocal; double NGlobal; double mu; double nu; //double fNoise;
    
    double dSfIf; double dSfIg; double dSgIg; double dSgIf; double dIfRf; double dIgRg;
    double dSfD; double dIfD; double dRfD; double dSgD; double dIgD; double dRgD; double dNfB; double dNgB; double Rf; double Rg;
    double newSf; double newSg; double newIf; double newIg; double newRf;  double newRg; double newNf; double newNg;
    
    double initSFocal; double initIFocal;
    
    @Override
    public void setInitParams()
    {
        
        //Set params to initial values and determine their order in params list
        betaFF = 0.48; params.add(betaFF); estParams.add(0);
        betaGF = 0.0 * betaFF; params.add(betaGF); fixedParams.add(1); //set betaGF and betaGF to zero for unstructured model
        betaGG = 0.24; params.add(betaGG); fixedParams.add(2);
        betaFG = betaGF; params.add(betaFG); fixedParams.add(3);
        
        alphaFocal = 0.068; params.add(alphaFocal); estParams.add(4);
        alphaGlobal = 0.0; params.add(alphaGlobal); fixedParams.add(5);
        deltaFocal = 0.455; params.add(deltaFocal); estParams.add(6);
        deltaGlobal = 0.455; params.add(deltaGlobal); fixedParams.add(7);
        NFocal = 10.0e06; params.add(NFocal); fixedParams.add(8);
        NGlobal = 25.0e06; params.add(NGlobal); fixedParams.add(9);
        mu = 1/(60.0*365.25); params.add(mu); fixedParams.add(10);
        nu = 1.0/7.0; params.add(nu); fixedParams.add(11);
        
        initSFocal = 0.2918; params.add(initSFocal); estParams.add(12); //estimated as the fraction, not the absolute number
        //initSGlobal = 0.5913; params.add(initSGlobal); fixedParams.add(13);
        
        //initIFocal = 0.00019436; params.add(initIFocal); estParams.add(13); // was 0.0005;
        //initIGlobal = 0.00028877;
                        
    }
    
    @Override
    public void setRandomGenerator() 
    {
      
        engine = new cern.jet.random.engine.MersenneTwister(new java.util.Date());
        poissDist = new Poisson(0.0, engine);
        normDist = new Normal(0.0, 1.0, engine);
        
    }
    
    @Override
    public void updateEstParams(DoubleArrayList newParams)
    {
        int p = estParams.size();
        int paramListIndex;
        for (int i = 0; i < p; i++) {
            paramListIndex = estParams.get(i);
            params.setQuick(paramListIndex, newParams.get(i));
        }
        
        //Update the estimated params
        betaFF = params.getQuick(0);
        //betaGF = params.getQuick(1);
        //betaGG = params.getQuick(2);
        //betaFG = betaGF; params.setQuick(3,betaGF);
        alphaFocal = params.getQuick(4);
        //alphaGlobal = alphaFocal; params.setQuick(5, alphaFocal);
        deltaFocal = params.getQuick(6);
        //deltaGlobal = deltaFocal; params.setQuick(7, deltaFocal);
        initSFocal = params.getQuick(12);
        //initIFocal = params.getQuick(13);
        
    }
    
    @Override
    public DoubleArrayList getInitialConditions() 
    {
        
        //Order init conditions 1) Sf 2) Sg 3) If 4) Ig 5) Nf 6) Ng
        
        currInits = new DoubleArrayList();
       
        //Compute equilibrium conditons without migration
        double R0Focal = betaFF / (mu + nu); //assuming no migration
        double R0Global = betaGG / (mu + nu); //assuming no migration
        
        double initS1 = NFocal/R0Focal; currInits.add(initS1);
        double initS2 = NGlobal/R0Global; currInits.add(initS2);
        double initI1 = mu * NFocal * (R0Focal - 1) / betaFF; currInits.add(initI1);
        double initI2 = mu * NGlobal * (R0Global - 1) / betaGG; currInits.add(initI2);
        currInits.add(NFocal);
        currInits.add(NGlobal);
        
        double dt = 0.25;
        double length = 365.25*100;
        currInits = this.fastSim(currInits, dt, length); //simulate to get endemic initial conditions for global pop
        
        //If estimating as a free parameter
        initS1 = initSFocal * NFocal; currInits.set(0, initS1); //estimating Sf init
        //initI1 = initIFocal * NFocal; currInits.set(2, initI1);
        
        //currInits.add(0.0); //dSfIf;
        //currInits.add(0.0); //dSfIg;
        //currInits.add(0.0); //dSgIg;
        //currInits.add(0.0); //dSgIf;
        //currInits.add(0.0); //IfRf;
        //currInits.add(0.0); //IgRg;
        //currInits.add(0.0); //NfB;
        //currInits.add(0.0); //NgB;        
        //currInits.add(0.0); //SfD;
        //currInits.add(0.0); //IfD;
        //currInits.add(0.0); //RfD;
        //currInits.add(0.0); //SgD;
        //currInits.add(0.0); //IgD;
        //currInits.add(0.0); //RgD;
        
        currInits.add(0.0); //cumulative incidence in focal pop
        currInits.add(0.0); //cumulative incidence in global pop
        
        return currInits;
        
    }
    
    @Override
    public double getPriorProb() 
    {
        double prior = 1.0;
        
        //Prior on betaFF - uniform between 0 and 10
        double betaFFPrior = 0.1;
        
        //All other priors - uniform between 0 and 1
        
        prior *= betaFFPrior;
        
        return Math.log(prior);
    }
    
    @Override
    public boolean checkParamConstraints()
    {
        boolean constraintCheckFail = false;
        
        int p = estParams.size();
        int paramListIndex;
        //Check if any params are negative
        for (int i = 0; i < p; ++i) {
            paramListIndex = estParams.get(i);
            if (params.get(paramListIndex) <= 0.0) {
                constraintCheckFail = true;
            }
        }
        
        if (params.get(6) > 1.0) {
            constraintCheckFail = true;
        }
        
        if (params.get(7) > 1.0) {
            constraintCheckFail = true;
        }
        
        if (params.get(12) > 1.0) {
            constraintCheckFail = true;
        }
        
        //if (params.get(13) > 1.0) {
            //constraintCheckFail = true;
        //}

        return constraintCheckFail;
    }
    
    @Override
    public DoubleMatrix3D updateStatesDeterministic(int jParticles, DoubleMatrix3D matrix, int xLocStart, DoubleArrayList xDtTimes, DoubleArrayList xSubTimes) 
    {
        
        //Spatial SIR with focal and global population
        //Order of state variables is 1) Sf 2) Sg 3) If 4) Ig 5) Nf 6) Ng
        
	for (int n = 0; n < xDtTimes.size(); ++n) {
	
            int xIndex = xLocStart + n;
            double dt = xDtTimes.getQuick(n);
            int nextIndex = xIndex + 1;
            
            double timeNow = xSubTimes.get(n);
            
            //Compute betaSeas for focal pop
            seasNowFocal = ((timeNow/365.25) + deltaFocal) - Math.floor(timeNow/365.25);
            seasHeightFocal = Math.cos(2*Math.PI * seasNowFocal);
            betaSeasFF = betaFF * (1 + alphaFocal*seasHeightFocal);
            betaSeasGF = betaGF * (1 + alphaFocal*seasHeightFocal); 
            
            //Compute betaSeas for global pop
            seasNowGlobal = ((timeNow/365.25) + deltaGlobal) - Math.floor(timeNow/365.25);
            seasHeightGlobal = Math.cos(2*Math.PI * seasNowGlobal);
            betaSeasGG = betaGG * (1 + alphaGlobal*seasHeightGlobal);
            betaSeasFG = betaFG * (1 + alphaGlobal*seasHeightGlobal);
            
            for (int j = 0; j < jParticles; ++j) {
                
                //double nScaler = 1.0;

                //INFECTIONS IN FOCAL POP FROM FOCAL POP
                dSfIf = (betaSeasFF * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(4,j,xIndex));
                
                //INFECTIONS IN FOCAL POP FROM GLOBAL POP
                dSfIg = (betaSeasGF * matrix.getQuick(0,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(4,j,xIndex));

                //INFECTIONS IN GLOBAL POP FROM GLOBAL POP
                dSgIg = (betaSeasGG * matrix.getQuick(1,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(5,j,xIndex));

                //INFECTIONS IN GLOBAL POP FROM FOCAL POP
                dSgIf = (betaSeasFG * matrix.getQuick(1,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(5,j,xIndex));
                
                //RECOVERIES
                dIfRf = (nu * matrix.getQuick(2,j,xIndex) * dt);
                dIgRg = (nu * matrix.getQuick(3,j,xIndex) * dt);
                
                //BIRTHS
                dNfB = (mu * matrix.getQuick(4,j,xIndex) * dt); //using actual birth rate
                dNgB = (mu * matrix.getQuick(5,j,xIndex) * dt); //usign actual birth rate
                
                //DEATHS IN FOCAL POP
                Rf = matrix.getQuick(4,j,xIndex) - matrix.getQuick(0,j,xIndex) - matrix.getQuick(2,j,xIndex);
                Rg = matrix.getQuick(5,j,xIndex) - matrix.getQuick(1,j,xIndex) - matrix.getQuick(3,j,xIndex);
                
                dSfD = (mu * matrix.getQuick(0,j,xIndex) * dt);
                dIfD = (mu * matrix.getQuick(2,j,xIndex) * dt);
                dRfD = (mu * Rf * dt);
                
                dSgD = (mu * matrix.getQuick(1,j,xIndex) * dt);
                dIgD = (mu * matrix.getQuick(3,j,xIndex) * dt);
                dRgD = (mu * Rg * dt);

                newSf = matrix.getQuick(0,j,xIndex) + dNfB - dSfIf - dSfIg - dSfD;
                newSg = matrix.getQuick(1,j,xIndex) + dNgB - dSgIg - dSgIf - dSgD;
                newIf = matrix.getQuick(2,j,xIndex) + dSfIf + dSfIg - dIfRf - dIfD;
                if (newIf < 0.0) {
                    newIf = 0.0;
                }
                newIg = matrix.getQuick(3,j,xIndex) + dSgIg + dSgIf - dIgRg - dIgD;
                if (newIg < 0.0) {
                    newIg = 0.0;
                }
                newRf = Rf + dIfRf - dRfD;
                newRg = Rg + dIgRg - dRgD;
                newNf = newSf + newIf + newRf;
                newNg = newSg + newIg + newRg;
                
                matrix.setQuick(0,j,nextIndex,newSf);
                matrix.setQuick(1,j,nextIndex,newSg);
                matrix.setQuick(2,j,nextIndex,newIf);
                matrix.setQuick(3,j,nextIndex,newIg);
                matrix.setQuick(4,j,nextIndex,newNf);
                matrix.setQuick(5,j,nextIndex,newNg);
                //matrix.setQuick(6,j,nextIndex,dSfIf);
                //matrix.setQuick(7,j,nextIndex,dSfIg);
                //matrix.setQuick(8,j,nextIndex,dSgIg);
                //matrix.setQuick(9,j,nextIndex,dSgIf);
                //matrix.setQuick(10,j,nextIndex,dIfRf);
                //matrix.setQuick(11,j,nextIndex,dIgRg);
                //matrix.setQuick(12,j,nextIndex,dNfB);
                //matrix.setQuick(13,j,nextIndex,dNgB);
                //matrix.setQuick(14,j,nextIndex,dSfD);
                //matrix.setQuick(15,j,nextIndex,dIfD);
                //matrix.setQuick(16,j,nextIndex,dRfD);
                //matrix.setQuick(17,j,nextIndex,dSgD);
                //matrix.setQuick(18,j,nextIndex,dIgD);
                //matrix.setQuick(19,j,nextIndex,dRgD);
                
                //Track cumulative incidence
                matrix.setQuick(6,j,nextIndex, (dSfIf + dSfIg));
                matrix.setQuick(7,j,nextIndex, (dSgIg + dSgIf));
                
            }	
        }
        return matrix;
    }
    
    public DoubleArrayList fastSim(DoubleArrayList inits, double dt, double endTime)
    {
        
        //Get sim times
        DoubleArrayList times = new DoubleArrayList();
        double timeNow = 0.0;
        times.add(timeNow);
        while (timeNow < endTime) {
            double timeNew = timeNow + dt;
            if (timeNew <= endTime) {
                times.add(timeNew);
                timeNow = timeNew;
            } else {
                times.add(endTime);
                timeNow = timeNew;
            }
        }
        
        //Set up array for state variables
        DoubleArrayList x = new DoubleArrayList();
        x.add(inits.get(0));
        x.add(inits.get(1));
        x.add(inits.get(2));
        x.add(inits.get(3));
        x.add(inits.get(4));
        x.add(inits.get(5));
        
	for (int n = 0; n < times.size(); ++n) {
            
            timeNow = times.getQuick(n);
            
            //Compute betaSeas for focal pop
            seasNowFocal = ((timeNow/365.25) + deltaFocal) - Math.floor(timeNow/365.25);
            seasHeightFocal = Math.cos(2*Math.PI * seasNowFocal);
            betaSeasFF = betaFF * (1 + alphaFocal*seasHeightFocal);
            betaSeasGF = betaGF * (1 + alphaFocal*seasHeightFocal); 
            
            //Compute betaSeas for global pop
            seasNowGlobal = ((timeNow/365.25) + deltaGlobal) - Math.floor(timeNow/365.25);
            seasHeightGlobal = Math.cos(2*Math.PI * seasNowGlobal);
            betaSeasGG = betaGG * (1 + alphaGlobal*seasHeightGlobal);
            betaSeasFG = betaFG * (1 + alphaGlobal*seasHeightGlobal);            
            
            //stndNormRand = stndNorm.nextDouble();
            //eta = stndNormRand / Math.sqrt(dt); //used for this particle and time only

            //INFECTIONS IN FOCAL POP FROM FOCAL POP
            //fTerm = fNoise * betaSeasFocal * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * eta * dt / matrix.getQuick(4, j, xIndex);
            dSfIf = (betaSeasFF * x.getQuick(0) * x.getQuick(2) * dt / x.getQuick(4)); // + fTerm;

            //INFECTIONS IN FOCAL POP FROM GLOBAL POP
            dSfIg = betaSeasGF * x.getQuick(0) * x.getQuick(3) * dt / x.getQuick(4);
 

            //INFECTIONS IN GLOBAL POP FROM GLOBAL POP
            //fTerm = fNoise * betaSeasGlobal * matrix.getQuick(1,j,xIndex) * matrix.getQuick(3,j,xIndex) * eta * dt / matrix.getQuick(5, j, xIndex);
            dSgIg = (betaSeasGG * x.getQuick(1) * x.getQuick(3) * dt / x.getQuick(5)); //+ fTerm;

            //INFECTIONS IN GLOBAL POP FROM FOCAL POP
            dSgIf = betaSeasFG * x.getQuick(1) * x.getQuick(2) * dt / x.getQuick(5);

            //RECOVERIES
            dIfRf = nu * x.getQuick(2) * dt;
            dIgRg = nu * x.getQuick(3) * dt;

            //BIRTHS
            dNfB = mu * x.getQuick(4) * dt;
            dNgB = mu * x.getQuick(5) * dt;

            //DEATHS IN FOCAL POP
            Rf = x.getQuick(4) - x.getQuick(0) - x.getQuick(2);
            Rg = x.getQuick(5) - x.getQuick(1) - x.getQuick(3);
            dSfD = mu * x.getQuick(0) * dt;
            dIfD = mu * x.getQuick(2) * dt;
            dRfD = mu * Rf * dt;

            dSgD = mu * x.getQuick(1) * dt;
            dIgD = mu * x.getQuick(3) * dt;
            dRgD = mu * Rg * dt;   

            newSf = x.getQuick(0) + dNfB - dSfIf - dSfIg - dSfD;
            newSg = x.getQuick(1) + dNgB - dSgIg - dSgIf - dSgD;
            newIf = x.getQuick(2) + dSfIf + dSfIg - dIfRf - dIfD;
            newIg = x.getQuick(3) + dSgIg + dSgIf - dIgRg - dIgD;
            newRf = Rf + dIfRf - dRfD;
            newRg = Rg + dIgRg - dRgD;
            newNf = newSf + newIf + newRf;
            newNg = newSg + newIg + newRg;

            x.setQuick(0,newSf);
            x.setQuick(1,newSg);
            x.setQuick(2,newIf);
            x.setQuick(3,newIg);
            x.setQuick(4,newNf);
            x.setQuick(5,newNg);
	
        }
        return x;
    }

       
}
