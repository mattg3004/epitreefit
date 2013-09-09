package epitreefit;

import cern.colt.list.*;
import cern.jet.random.*;
import cern.colt.matrix.*;

/**
 * This is the model used in MBE paper for vectored and combined vector/spatial model
 * @author David
 */

public class EpiModelDengueVectored extends EpiModel {

    Poisson poissDist;
    
    double seasNowFocal; double seasNowGlobal; double seasHeightFocal; double seasHeightGlobal; double betaSeasFF; double betaSeasGF; double betaSeasGG; double betaSeasFG;
    double betaHV; double betaVH; double betaGF; double betaGG; double betaFG; double alphaFocal; double alphaGlobal; double deltaFocal; double deltaGlobal;
    double NFocal; double NGlobal; double mu; double nu; //double fNoise;
    
    double dSfIf; double dSfIg; double dSgIg; double dSgIf; double dIfRf; double dIgRg;
    double dSfD; double dIfD; double dRfD; double dSgD; double dIgD; double dRgD; double dNfB; double dNgB; double Rf; double Rg;
    double newSf; double newSg; double newIf; double newIg; double newRf;  double newRg; double newNf; double newNg;
    
    double M; double muVec;
    double B; double seasB; double seasM;
    
    double initSFocal; double initIFocal;
    
    @Override
    public void setInitParams()
    {
        //Set params to initial values and determine their order in params list
        //betaFF = 0.22; params.add(betaFF); estParams.add(0);
        betaHV = 0.2770; params.add(betaHV); estParams.add(0); //transmission rate from humans to mosquitoes in focal pop
        betaVH = betaHV; params.add(betaVH); fixedParams.add(1);
        betaGF = 0.0; params.add(betaGF); fixedParams.add(2); //set betaGF and betaFG to zero for model without spatial structure
        betaGG = 0.1595; params.add(betaGG); fixedParams.add(3);
        betaFG = betaGF; params.add(betaFG); fixedParams.add(4);
        
        alphaFocal = 0.11; params.add(alphaFocal); estParams.add(5);
        alphaGlobal = 0.0; params.add(alphaGlobal); fixedParams.add(6);
        deltaFocal = 6.0/12.0; params.add(deltaFocal); estParams.add(7);
        deltaGlobal = 6.0/12.0; params.add(deltaGlobal); fixedParams.add(8);
        NFocal = 10.0e06; params.add(NFocal); fixedParams.add(9);
        NGlobal = 25.0e06; params.add(NGlobal); fixedParams.add(10);
        mu = 1/(60*365.25); params.add(mu); fixedParams.add(11);
        nu = 1.0/7.0; params.add(nu); fixedParams.add(12);
        
        //Additional vector params;
        M = 0.65; params.add(M); estParams.add(13);
        muVec = 1.0 / 7.0; params.add(muVec); fixedParams.add(14);
        
        //Init conditions
        initSFocal = 0.4041; params.add(initSFocal); estParams.add(15); //estimated as the fraction, not the absolute number
        initIFocal = 0.0005; params.add(initIFocal); fixedParams.add(16); //estimated as the fraction, not the absolute number
                        
    }
    
    @Override
    public void setRandomGenerator() 
    {
      
        cern.jet.random.engine.RandomEngine engine = new cern.jet.random.engine.MersenneTwister(new java.util.Date());
        poissDist = new Poisson(0.0, engine);
        
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
        betaHV = params.getQuick(0);
        betaVH = params.getQuick(0); params.setQuick(1, betaVH);
        //betaGF = params.getQuick(2);
        //betaGG = params.getQuick(3);
        //betaFG = params.getQuick(2); params.setQuick(4, betaFG);
        alphaFocal = params.getQuick(5);
        deltaFocal = params.getQuick(7);
        
        M = params.getQuick(13);
        initSFocal = params.getQuick(15);
        initIFocal = params.getQuick(16);
        
    }
    
    @Override
    public DoubleArrayList getInitialConditions() 
    {
        
        //Order init conditions 1) Sf 2) Sg 3) If 4) Ig 5) Nf 6) Ng
        
        DoubleArrayList currInits = new DoubleArrayList();
       
        //Compute equilibrium conditons without migration
        double R0Focal = (betaHV*betaVH * M) / (muVec * (mu + nu));
        double betaFocal = R0Focal * (mu + nu); 
        double R0Global = betaGG / (mu + nu); //assuming no migration
        
        double initS1 = NFocal/R0Focal;
        double initS2 = NGlobal/R0Global;
        double initI1 = mu * NFocal * (R0Focal - 1) / (betaFocal); 
        double initI2 = mu * NGlobal * (R0Global - 1) / betaGG;
        
        double initSv = (muVec * M * NFocal) / ((betaHV/NFocal) * initI1 + muVec);
        double initIv = (betaHV * initSv * (initI1/NFocal)) / muVec;
        
        currInits.add(initS1);
        currInits.add(initS2);
        currInits.add(initSv);
        currInits.add(initI1);
        currInits.add(initI2);
        currInits.add(initIv);
        currInits.add(NFocal);
        currInits.add(NGlobal);
        
        double dt = 0.25;
        double length = 365.25*200;
        currInits = this.fastSim(currInits, dt, length);
        
        //If estimating init S1
        initS1 = initSFocal * NFocal; currInits.set(0, initS1); //estimating Sf init
        
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
        
        //Need to add trackers for If and Ig (or Iv)
        currInits.add(0.0); //cumulative incidence in focal pop
        currInits.add(0.0); //cumulative incidence in global pop
        
        return currInits;
        
    }
    
    @Override
    public double getPriorProb() 
    {
        double prior = 1.0;
        
        //Prior on betaVH and betaHV - uniform between 0 and 10
        double betaVHPrior = 0.1;
        
        //Prior on M - uniform between 0 and 100
        double mPrior = 0.01;
        
        //All other priors - uniform between 0 and 1
        
        prior *= betaVHPrior;
        prior *= mPrior;
        
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
        
        //If estimating init S1
        if (params.get(15) > 1.0) {
            constraintCheckFail = true;
        }
        
        //deltaFocal has to be less than 365.25
        if (params.get(7) > (1.0)) { //this means delta is less than one year
            constraintCheckFail = true;
        }
        
        
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
            //betaSeasFF = betaFF * (1 + alphaFocal*seasHeightFocal);
            betaSeasGF = betaGF * (1 + alphaFocal*seasHeightFocal); 
            
            //If estimating seasonal M
            //B = muVec * M * NFocal;
            seasM = M * (1 + alphaFocal*seasHeightFocal);
            
            //Compute betaSeas for global pop
            seasNowGlobal = ((timeNow/365.25) + deltaGlobal) - Math.floor(timeNow/365.25);
            seasHeightGlobal = Math.cos(2*Math.PI * seasNowGlobal);
            betaSeasGG = betaGG * (1 + alphaGlobal*seasHeightGlobal);
            betaSeasFG = betaFG * (1 + alphaGlobal*seasHeightGlobal);
            
            for (int j = 0; j < jParticles; ++j) {
                
                //Get quasi-equilibrium vector states
                //double vecI = seasM * NFocal * (betaFF/NFocal) * muVec * matrix.getQuick(2,j,xIndex) / ((betaFF/NFocal)*matrix.getQuick(2,j,xIndex)*muVec + muVec*muVec);
                //double vecI = (seasB * betaHV * (matrix.getQuick(2,j,xIndex) / NFocal)) / ((betaHV/NFocal)*matrix.getQuick(2,j,xIndex)*muVec + muVec*muVec);

                //INFECTIONS IN FOCAL HUMANS POP FROM FOCAL VECTORS POP
                //dSfIf = (betaSeasFF * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(4,j,xIndex));
                dSfIf = (betaVH * matrix.getQuick(0,j,xIndex) * matrix.getQuick(5,j,xIndex) * dt / matrix.getQuick(6,j,xIndex));
                
                //INFECTIONS IN FOCAL VECTORS FROM FOCAL HUMANS
                double dSvIf = (betaHV * matrix.getQuick(2,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(6,j,xIndex));
                
                //INFECTIONS IN FOCAL POP FROM GLOBAL POP
                dSfIg = (betaSeasGF * matrix.getQuick(0,j,xIndex) * matrix.getQuick(4,j,xIndex) * dt / matrix.getQuick(6,j,xIndex));

                //INFECTIONS IN GLOBAL POP FROM GLOBAL POP
                dSgIg = (betaSeasGG * matrix.getQuick(1,j,xIndex) * matrix.getQuick(4,j,xIndex) * dt / matrix.getQuick(7,j,xIndex));

                //INFECTIONS IN GLOBAL POP FROM FOCAL POP
                dSgIf = (betaSeasFG * matrix.getQuick(1,j,xIndex) * matrix.getQuick(3,j,xIndex) * dt / matrix.getQuick(7,j,xIndex));
                
                //RECOVERIES
                dIfRf = (nu * matrix.getQuick(3,j,xIndex) * dt);
                dIgRg = (nu * matrix.getQuick(4,j,xIndex) * dt);
                
                //BIRTHS
                dNfB = (mu * matrix.getQuick(6,j,xIndex) * dt);
                dNgB = (mu * matrix.getQuick(7,j,xIndex) * dt);
                
                //DEATHS IN HUMANS
                Rf = matrix.getQuick(6,j,xIndex) - matrix.getQuick(0,j,xIndex) - matrix.getQuick(3,j,xIndex);
                Rg = matrix.getQuick(7,j,xIndex) - matrix.getQuick(1,j,xIndex) - matrix.getQuick(4,j,xIndex);
                
                dSfD = (mu * matrix.getQuick(0,j,xIndex) * dt);
                dIfD = (mu * matrix.getQuick(3,j,xIndex) * dt);
                dRfD = (mu * Rf * dt);
                
                dSgD = (mu * matrix.getQuick(1,j,xIndex) * dt);
                dIgD = (mu * matrix.getQuick(4,j,xIndex) * dt);
                dRgD = (mu * Rg * dt);
                
                //DEATHS IN MOSQUITOES
                double dIvD = muVec * matrix.getQuick(5,j,xIndex) * dt;
                
                //UPDATES
                newSf = matrix.getQuick(0,j,xIndex) + dNfB - dSfIf - dSfIg - dSfD;
                newSg = matrix.getQuick(1,j,xIndex) + dNgB - dSgIg - dSgIf - dSgD;
                newIf = matrix.getQuick(3,j,xIndex) + dSfIf + dSfIg - dIfRf - dIfD;
                if (newIf < 0.0) {
                    newIf = 0.0;
                }
                newIg = matrix.getQuick(4,j,xIndex) + dSgIg + dSgIf - dIgRg - dIgD;
                if (newIg < 0.0) {
                    newIg = 0.0;
                }
                newRf = Rf + dIfRf - dRfD;
                newRg = Rg + dIgRg - dRgD;
                newNf = newSf + newIf + newRf;
                newNg = newSg + newIg + newRg;
                
                double newIv = matrix.getQuick(5,j,xIndex) + dSvIf - dIvD;
                double newNv = seasM * matrix.getQuick(6,j,xIndex);
                double newSv = newNv - newIv;
                
                matrix.setQuick(0,j,nextIndex,newSf);
                matrix.setQuick(1,j,nextIndex,newSg);
                matrix.setQuick(2,j,nextIndex,newSv);
                matrix.setQuick(3,j,nextIndex,newIf);
                matrix.setQuick(4,j,nextIndex,newIg);
                matrix.setQuick(5,j,nextIndex,newIv);
                matrix.setQuick(6,j,nextIndex,newNf);
                matrix.setQuick(7,j,nextIndex,newNg);
                //matrix.setQuick(8,j,nextIndex,dSfIf);
                //matrix.setQuick(9,j,nextIndex,dSfIg);
                //matrix.setQuick(10,j,nextIndex,dSgIg);
                //matrix.setQuick(11,j,nextIndex,dSgIf);
                //matrix.setQuick(12,j,nextIndex,dIfRf);
                //matrix.setQuick(13,j,nextIndex,dIgRg);
                //matrix.setQuick(14,j,nextIndex,dNfB);
                //matrix.setQuick(15,j,nextIndex,dNgB);
                //matrix.setQuick(16,j,nextIndex,dSfD);
                //matrix.setQuick(17,j,nextIndex,dIfD);
                //matrix.setQuick(18,j,nextIndex,dRfD);
                //matrix.setQuick(19,j,nextIndex,dSgD);
                //matrix.setQuick(20,j,nextIndex,dIgD);
                //matrix.setQuick(21,j,nextIndex,dRgD);
                
                //Track cumulative incidence
                matrix.setQuick(8,j,nextIndex, (dSfIf + dSfIg));
                matrix.setQuick(9,j,nextIndex, (dSgIg + dSgIf));
                
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
        x.add(inits.get(6));
        x.add(inits.get(7));
        
	for (int n = 0; n < times.size(); ++n) {
            
            timeNow = times.getQuick(n);
            
            //Compute betaSeas for focal pop
            seasNowFocal = ((timeNow/365.25) + deltaFocal) - Math.floor(timeNow/365.25);
            seasHeightFocal = Math.cos(2*Math.PI * seasNowFocal);
            //betaSeasFF = betaFF * (1 + alphaFocal*seasHeightFocal);
            betaSeasGF = betaGF * (1 + alphaFocal*seasHeightFocal);
            
            //If estimating seasonal M
            //B = muVec * M * NFocal;
            //seasB = B * (1 + alphaFocal*seasHeightFocal);
            seasM = M * (1 + alphaFocal*seasHeightFocal);
            
            //Compute betaSeas for global pop
            seasNowGlobal = ((timeNow/365.25) + deltaGlobal) - Math.floor(timeNow/365.25);
            seasHeightGlobal = Math.cos(2*Math.PI * seasNowGlobal);
            betaSeasGG = betaGG * (1 + alphaGlobal*seasHeightGlobal);
            betaSeasFG = betaFG * (1 + alphaGlobal*seasHeightGlobal);            
            
            //INFECTIONS IN FOCAL HUMANS POP FROM FOCAL VECTORS POP
            //dSfIf = (betaSeasFF * matrix.getQuick(0,j,xIndex) * matrix.getQuick(2,j,xIndex) * dt / matrix.getQuick(4,j,xIndex));
            dSfIf = (betaVH * x.getQuick(0) * x.getQuick(5) * dt / x.getQuick(6));

            //INFECTIONS IN FOCAL VECTORS FROM FOCAL HUMANS
            double dSvIf = (betaHV * x.getQuick(2) * x.getQuick(3) * dt / x.getQuick(6));

            //INFECTIONS IN FOCAL POP FROM GLOBAL POP
            dSfIg = (betaSeasGF * x.getQuick(0) * x.getQuick(4) * dt / x.getQuick(6));

            //INFECTIONS IN GLOBAL POP FROM GLOBAL POP
            dSgIg = (betaSeasGG * x.getQuick(1) * x.getQuick(4) * dt / x.getQuick(7));

            //INFECTIONS IN GLOBAL POP FROM FOCAL POP
            dSgIf = (betaSeasFG * x.getQuick(1) * x.getQuick(3) * dt / x.getQuick(7));

            //RECOVERIES
            dIfRf = (nu * x.getQuick(3) * dt);
            dIgRg = (nu * x.getQuick(4) * dt);

            //BIRTHS
            dNfB = (mu * x.getQuick(6) * dt);
            dNgB = (mu * x.getQuick(7) * dt);

            //DEATHS IN HUMANS
            Rf = x.getQuick(6) - x.getQuick(0) - x.getQuick(3);
            Rg = x.getQuick(7) - x.getQuick(1) - x.getQuick(4);

            dSfD = (mu * x.getQuick(0) * dt);
            dIfD = (mu * x.getQuick(3) * dt);
            dRfD = (mu * Rf * dt);

            dSgD = (mu * x.getQuick(1) * dt);
            dIgD = (mu * x.getQuick(4) * dt);
            dRgD = (mu * Rg * dt);

            //DEATHS IN MOSQUITOES
            double dIvD = muVec * x.getQuick(5) * dt;

            //UPDATES
            newSf = x.getQuick(0) + dNfB - dSfIf - dSfIg - dSfD;
            newSg = x.getQuick(1) + dNgB - dSgIg - dSgIf - dSgD;
            newIf = x.getQuick(3) + dSfIf + dSfIg - dIfRf - dIfD;
            if (newIf < 0.0) {
                newIf = 0.0;
            }
            newIg = x.getQuick(4) + dSgIg + dSgIf - dIgRg - dIgD;
            if (newIg < 0.0) {
                newIg = 0.0;
            }
            newRf = Rf + dIfRf - dRfD;
            newRg = Rg + dIgRg - dRgD;
            newNf = newSf + newIf + newRf;
            newNg = newSg + newIg + newRg;

            double newIv = x.getQuick(5) + dSvIf - dIvD;
            double newNv = seasM * x.getQuick(6);
            double newSv = newNv - newIv;

            x.setQuick(0,newSf);
            x.setQuick(1,newSg);
            x.setQuick(2,newSv);
            x.setQuick(3,newIf);
            x.setQuick(4,newIg);
            x.setQuick(5,newIv);
            x.setQuick(6,newNf);
            x.setQuick(7,newNg);
	
        }
        return x;
    }
       
}
