package epitreefit;

import cern.colt.list.*;
import cern.colt.matrix.*;
//import cern.jet.random.*;
//import java.util.ArrayList;

/**
 * Simulate state trajectories and compute coalescent likelihoods
 * @author David
 */
public class ParticleFilter
{
    int jParticles; //number of particles
    double logLikelihood;
    StateTrajectory x;
    XTrajectory xTrajSample;
    LineageStateProbs stateProbs;
    DoubleFactory2D factory2D = DoubleFactory2D.dense;
    DoubleFactory3D factory3D = DoubleFactory3D.dense;
    
    public void runFilter(int jNum, EpiModel epi, ZVectors dataZ, TreeNode[] tree, StructCoalModel coal, double startTime, double endTime)
    {
	
        //Forward simulate the state trajectories from the EpiModel and then compute coalescent likelihood while integrating over lineage states
        
        jParticles = jNum;
        
        int zLocStart = dataZ.absoluteTimes.indexOf(startTime);
        int zLocEnd = dataZ.absoluteTimes.indexOf(endTime);
        xTrajSample = new XTrajectory(); //holds sampled trajectory obtained from particle filter
        xTrajSample.zLocStart = zLocStart;
        xTrajSample.zLocEnd = zLocEnd;
        DoubleArrayList xTimes = (DoubleArrayList) dataZ.absoluteTimes.partFromTo(zLocStart, zLocEnd);
	int totalTimes = xTimes.size();
        
        DoubleArrayList filterTimes = new DoubleArrayList();
        filterTimes.add(startTime);
        filterTimes.add(endTime);
        int filterTotalTimes = filterTimes.size();
        
	//Set up x matrix
        x = new StateTrajectory();
        x.getMatrix(jParticles, totalTimes, epi.getInitialConditions(), factory3D);
        
	//Set up w matrix
        Weights w = new Weights();
        w.getMatrix(jParticles, totalTimes, factory2D);
        
        //Run simulation from EpiModel
        double timeStart; double timeEnd; 	
        int xLocStart = 0; int xLocEnd = 0;
	for (int time = 0; time < (filterTotalTimes - 1); ++time) {
            
            timeStart = filterTimes.get(time); //absolute start time
            timeEnd = filterTimes.get(time+1); //absolute end time
            xLocStart = xTimes.indexOf(timeStart); //index of where to start in xTimes
            xLocEnd = xTimes.indexOf(timeEnd); //index of where to end in xTimes
            DoubleArrayList xSubTimes = (DoubleArrayList) xTimes.partFromTo(xLocStart, xLocEnd);
            DoubleArrayList xDtTimes = new DoubleArrayList();
            for (int i = 1; i < xSubTimes.size(); ++i) {
                double xDiff = xSubTimes.get(i) - xSubTimes.get(i-1);
		xDtTimes.add(xDiff);
            }

            //Update particle states
            x.propogateParticlesDeterministic(epi, xLocStart, xDtTimes, xSubTimes);
	}
        
        stateProbs = new LineageStateProbs();
        stateProbs.getMatrix(jParticles, tree.length, coal.Y.size(), factory3D);
        
        //Compute lineage state probabilities and coalescent likelihoods
        int backCounter = 0;
        for (int time = (totalTimes-1); time > 0; --time) {
            
            int zLocEndInterval = zLocEnd - backCounter; //zLoc at end of interval (time = t+1)
            int zLocStartInterval = zLocEndInterval - 1;
            double timeEndInterval = dataZ.absoluteTimes.get(zLocEndInterval); 
            double timeStartInterval = dataZ.absoluteTimes.get(zLocStartInterval);
            
            //Get new samples (if samples)
            if (dataZ.omegaEvents.get(zLocEndInterval).contains(2)) { //if this interval begins with a sampling event
                stateProbs.addSamples(dataZ, zLocEndInterval, tree);
            }
            
            //Update weights
            double dtTime = dataZ.absoluteTimes.get(zLocEndInterval) - dataZ.absoluteTimes.get(zLocStartInterval);
            w.updateWeightsNoCoal(x.matrix.viewColumn(time), stateProbs, epi.params, time, coal, dtTime, timeEndInterval);
            
            //Update lineageStateProbs
            stateProbs.updateProbs(x.matrix.viewColumn(time), epi.params, time, coal, dtTime, timeEndInterval);
            
            //Update weights (and lineage state probs) after coalescent events
            if (dataZ.omegaEvents.get(zLocStartInterval).contains(1)) {
                w.updateWeightsCoal(dataZ, zLocStartInterval, tree, x.matrix.viewColumn(time-1), stateProbs, epi.params, time, coal, timeStartInterval);
            }

            backCounter++;
        }
        
	logLikelihood = w.computeLogLikelihood();
        if (logLikelihood == Double.POSITIVE_INFINITY) { //this should never happen but...
            logLikelihood = Double.NEGATIVE_INFINITY;
            System.out.println("Log likelihood evaluated as INFINITY");
        }
        
        xTrajSample.full = x.matrix.viewRow(0); //get xTraj sample
        
    }//End Method

}
