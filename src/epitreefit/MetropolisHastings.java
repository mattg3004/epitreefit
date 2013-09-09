package epitreefit;

import cern.jet.random.*;
import cern.colt.list.*;
import cern.colt.matrix.*;
import java.util.ArrayList;
import java.io.*; //IOException;

/**
 * Main MH sampling algorithm
 * @author David
 */

public class MetropolisHastings
{
    public void runMCMC(int iterations, int particles, int outputFreq, TreeNode[] tree, double startTime, double endTime, double dt, ProposalDensity thetaProposal, EpiModel epi, StructCoalModel coal, String fileNameStem, String fileNamePostFix, boolean verbose) throws IOException
    {

        //Theta sample output file for model parameters
        SampleWriter sampleWriter = new SampleWriter();
        String paramFileName = fileNameStem + "_params_" + fileNamePostFix;
        FileWriter thetaFileWriter = new FileWriter(paramFileName);
        PrintWriter thetaPrintWriter = new PrintWriter(thetaFileWriter);
        
        //xTraj sample output file for state variables
        String xTrajFileName = fileNameStem + "_xTrajs_" + fileNamePostFix;
        FileWriter xTrajFileWriter = new FileWriter(xTrajFileName);
        PrintWriter xTrajPrintWriter = new PrintWriter(xTrajFileWriter);

        //Likelihood output file
        String likeFileName = fileNameStem + "_likes_" + fileNamePostFix;
        FileWriter likeFileWriter = new FileWriter(likeFileName);
        PrintWriter likePrintWriter = new PrintWriter(likeFileWriter);
        
        //Write init params to file
        DoubleArrayList thetaEstNow = epi.getEstParams(); //DONT NEED THIS
        DoubleArrayList thetaAllNow = thetaEstNow.copy();
        DoubleArrayList thetaFixed = epi.getFixedParams();
        thetaAllNow.addAllOf(thetaFixed);
        sampleWriter.samplesToString(thetaEstNow);
        thetaPrintWriter.print(sampleWriter.sampleString); thetaPrintWriter.flush(); thetaPrintWriter.println();

        //Set ZVectors - this is the main data structure used to compute the coalescent likelihood
        ZVectors dataZNow = new ZVectors(); //Vector data structure for tree and observation data
        dataZNow.setZVectors(tree, startTime, endTime, dt);
        startTime = dataZNow.startTime;
        endTime = dataZNow.endTime;
        
        ParticleFilter filter = new ParticleFilter(); //computes model likelihoods (this is a holdover from the particle filtering algorithms in Rasmussen et al. 2011)
        filter.runFilter(particles, epi, dataZNow, tree, coal, startTime, endTime);
        XTrajectory xTrajNow = filter.xTrajSample; //this should not be updated after this point!!
        double pGNow = filter.logLikelihood; //retrieve log likelihood
        double pThetaNow = epi.getPriorProb(); //get (log) prior prob of params in \theta

        //Initialize xTrajSamples for xTraj samples
        XTrajSamples xTrajs = new XTrajSamples();
        double sampleFreq = 30.0; //do this in weeks to align with sim times
        int sampleVar1 = 6; //incidence in HCMC
        int sampleVar2 = 7; //incidence in non-HCMC
        
        //If sampling prevalence from xTrajs
//        xTrajs.getSampleTimes(startTime, endTime, sampleFreq);
//        DoubleArrayList xTrajSampleNow = xTrajs.getSample(sampleVar1, dataZNow, xTrajNow);
//        sampleWriter.samplesToString(xTrajs.sampleTimes);
//        xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();
//        sampleWriter.samplesToString(xTrajSampleNow);
//        xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();
//        xTrajSampleNow = xTrajs.getSample(sampleVar2, dataZNow, xTrajNow);
//        sampleWriter.samplesToString(xTrajSampleNow);
//        xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();

        //If sampling incidence counts from xTrajs
        xTrajs.getSampleTimes(startTime, endTime, sampleFreq);
        sampleWriter.samplesToString(xTrajs.sampleTimes);
        xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();

        DoubleArrayList xTrajSampleNow = xTrajs.getSampleIncidence(sampleVar1, dataZNow, xTrajNow);
        sampleWriter.samplesToString(xTrajSampleNow);
        xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();

        xTrajSampleNow = xTrajs.getSampleIncidence(sampleVar2, dataZNow, xTrajNow);
        sampleWriter.samplesToString(xTrajSampleNow);
        xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();

        //Write out starting likelihoods
        DoubleArrayList likeSamples = new DoubleArrayList();
        likeSamples.add(pGNow); likeSamples.add(pThetaNow);
        sampleWriter.samplesToString(likeSamples);
        likePrintWriter.print(0 + ", "); likePrintWriter.print(sampleWriter.sampleString); likePrintWriter.flush(); likePrintWriter.println();

        //Set up autoOptimization routines for proposal density
        AutoOptimizer optimizer = new AutoOptimizer();
        optimizer.autoOptimize = true; //will optimze proposal covariance structure if true
        if (optimizer.autoOptimize) {
            ArrayList<Integer> runPoints = new ArrayList<Integer>();
            runPoints.add(200); runPoints.add(500); runPoints.add(1000); runPoints.add(2000);
            optimizer.setRunPoints(runPoints, thetaEstNow.size());
            optimizer.getSampleMatrix();
        }

        //Initialize MCMC timer and counters
        double startClockTime = System.currentTimeMillis();
        int thetaProposalCount = 0;
        int thetaAcceptCount = 0;
        int thinCounter = 0;
        
        //Run MCMC
        for (int n = 1; n <= iterations; ++n) {

                thinCounter++;
                if (verbose) {
                    System.out.println("MCMC iteration = " + n);
                }
                    
                //Propose and accept/reject new \theta* with MH step
                thetaProposalCount++;
                DoubleArrayList thetaEstNew = thetaProposal.nextProposal(thetaEstNow, epi);
                DoubleArrayList thetaAllNew = thetaEstNew.copy();
                thetaAllNew.addAllOf(thetaFixed);
                epi.updateEstParams(thetaEstNew); //This isn't actually necessary since parameters get updated automatically

                filter.runFilter(particles, epi, dataZNow, tree, coal, startTime, endTime);
                double pGNew = filter.logLikelihood;
                double pThetaNew = epi.getPriorProb();

                if (verbose) {
                    System.out.println("thetaNow: " + thetaEstNow);
                    System.out.println("thetaNew " + thetaEstNew);
                    System.out.println("pGNow: " + pGNow);
                    System.out.println("pGNew: " + pGNew);
                }

                double probAcceptTheta = Math.exp((pGNew + pThetaNew) - (pGNow + pThetaNow));
                if (probAcceptTheta >= Math.random()) {
                    thetaAllNow = thetaAllNew.copy();
                    thetaEstNow = thetaEstNew.copy();
                    pGNow = pGNew;
                    pThetaNow = pThetaNew;
                    xTrajNow = filter.xTrajSample;
                    if (verbose) {
                        System.out.println("Accepted theta proposal");
                    }
                    thetaAcceptCount++;

                } else {
                    epi.updateEstParams(thetaEstNow);
                }

                if (thinCounter == outputFreq) {

                    sampleWriter.samplesToString(thetaEstNow);
                    thetaPrintWriter.print(sampleWriter.sampleString); thetaPrintWriter.flush(); thetaPrintWriter.println();

                    likeSamples = new DoubleArrayList();
                    likeSamples.add(pGNow); likeSamples.add(pThetaNow);
                    sampleWriter.samplesToString(likeSamples);
                    likePrintWriter.print(n + ", "); likePrintWriter.print(sampleWriter.sampleString); likePrintWriter.flush(); likePrintWriter.println();
                    
                    //For incidence 
                    xTrajSampleNow = xTrajs.getSampleIncidence(sampleVar1, dataZNow, xTrajNow);
                    sampleWriter.samplesToString(xTrajSampleNow);
                    xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();
                    xTrajSampleNow = xTrajs.getSampleIncidence(sampleVar2, dataZNow, xTrajNow);
                    sampleWriter.samplesToString(xTrajSampleNow);
                    xTrajPrintWriter.print(sampleWriter.sampleString); xTrajPrintWriter.flush(); xTrajPrintWriter.println();
                    
                    thinCounter = 0;
                }

                if (optimizer.autoOptimize) {   
                    optimizer.addSample(thetaEstNow);
                    if (optimizer.runPoints.contains(n)) {
                        DoubleMatrix2D newCovMatrix = optimizer.autoOptimize(); 
                        thetaProposal.covMatrix = newCovMatrix;
                    }
                }
                
                if (verbose) {
                    System.out.println();
                }

        }
        //MCMC statistics
        double elapsedClockTime = System.currentTimeMillis() - startClockTime;
        double elapsedTimeSecs = elapsedClockTime / 1000;
        System.out.println("Total run time was:" + elapsedTimeSecs + "seconds");
        double iters = (double) iterations;
        double avgTimePerIteration = elapsedTimeSecs/ iters;
        System.out.println("Average time per iteration: " + avgTimePerIteration + " secs");

        System.out.println("Theta acceptance rate:");
        double accepts = (double) thetaAcceptCount;
        double its = (double) thetaProposalCount;
        double acceptRate = accepts / its;
        System.out.println(acceptRate);

        thetaPrintWriter.close();
        thetaFileWriter.close();
        likePrintWriter.close();
        likeFileWriter.close();
        xTrajPrintWriter.close();
        xTrajFileWriter.close();

    }//END method
    
}
