package epitreefit;

import cern.colt.list.*;
import cern.colt.matrix.*;
import java.util.ArrayList;

/**
 *
 * @author David
 */
public class LineageStateProbs {
    
    int particles;
    int states;
    int lineages;
    DoubleMatrix3D matrix;
    ArrayList<Integer> activeLineages = new ArrayList<Integer>();
    //For debugging:
    ArrayList<Integer> lineagesRemoved = new ArrayList<Integer>();
    ArrayList<Integer> lineagesAdded = new ArrayList<Integer>();
    
    public void getMatrix(int jParticles, int lines, int sts, DoubleFactory3D factory3D) {
        
        states = sts;
        particles = jParticles;
        lineages = lines;
        matrix = factory3D.make(particles, lineages, sts);
        
    }
    
    public void addSamples(ZVectors dataZ, int zLoc, TreeNode[] tree) 
    {
        
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        //double eventTime = dataZ.absoluteTimes.get(zLoc);
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event) == 2) { //if sampling event
                
                int newLineage = dataZ.nodePointers.get(zLoc).get(event); //indexed from zero
                activeLineages.add(newLineage);
                lineagesAdded.add(newLineage);
                int newState = tree[newLineage].nodeLabel - 1; //nodeLabels are indexed from one
                for (int j = 0; j < particles; j++) {
                    if (newState == -1) {
                        System.out.println("Sampled lineage not in defined state");
                    }
                    matrix.set(j, newLineage, newState, 1.0); //set new state prob to 1.0
                }
            }
        }
        
    }
    
    //Not used for dengue analysis
    public void addSamplesWithPriors(ZVectors dataZ, int zLoc, TreeNode[] tree) 
    {
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        double statePrior = 0.0;
        //double eventTime = dataZ.absoluteTimes.get(zLoc);
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event) == 2) { //if sampling event
                
                int newLineage = dataZ.nodePointers.get(zLoc).get(event); //indexed from zero
                activeLineages.add(newLineage);
                lineagesAdded.add(newLineage);
                for (int j = 0; j < particles; j++) {
                    for (int s = 0; s < states; s++) {
                        statePrior = tree[newLineage].statePriors.get(s);
                        matrix.setQuick(j, newLineage, s, statePrior); //set new state prob to 1.0
                    }
                }
            }
        }
        
    }
    
    //Update lineage state probabilities
    public void updateProbs(DoubleMatrix2D xCurr, DoubleArrayList theta, int time, StructCoalModel coal, double dtTime, double currTime) 
    {
        
        int numLineages = activeLineages.size();
        
        //For each particle
        for (int x = 0; x < particles; x++) { //use x's for particle indexes to avoid confusion with lineage indexes i and j
        
            coal.updateF(xCurr.viewColumn(x), theta, currTime);
            coal.updateG(xCurr.viewColumn(x), theta);
            coal.updateY(xCurr.viewColumn(x));
            
            //Precompute A(l) terms by summing probs in matrix
            DoubleArrayList lineagesInAByState = this.getLineagesByState(x);
            
            double omegaIn = 0.0; double omegaOut = 0.0;
            double omegaGain; double omegaLoss;
            
            //double lambdaSum = 0.0; //total rate of coalescence across all pairs of lineages
            for (int lin = 0; lin < numLineages; lin++) { //for all lineages i
                
                int i = activeLineages.get(lin);
                DoubleArrayList newLineageProbs = new DoubleArrayList();
                
                for (int k = 0; k < states; k++) {
                    
                    omegaGain = 0.0;
                    omegaLoss = 0.0;
                    
                    //Solve the master equation as would an ODE with discrete time step
                    for (int l = 0; l < states; l++) {    
                        if (l != k) {
                            omegaIn = ((coal.G.getQuick(k, l)/coal.Y.getQuick(l)) + (coal.F.getQuick(k, l)/coal.Y.getQuick(l))*((coal.Y.getQuick(k) - lineagesInAByState.getQuick(k))/coal.Y.getQuick(k))) * dtTime * matrix.get(x,i,l);
                            if (!Double.isNaN(omegaIn)) {
                                omegaGain += omegaIn;
                            }
                            omegaOut = ((coal.G.getQuick(l, k)/coal.Y.getQuick(k)) + (coal.F.getQuick(l, k)/coal.Y.getQuick(k))*((coal.Y.getQuick(l) - lineagesInAByState.getQuick(l))/coal.Y.getQuick(l))) * dtTime * matrix.get(x,i,k);
                            if (!Double.isNaN(omegaOut)) {
                                omegaLoss += omegaOut;
                            }
                        }
                    }
                    
                    double newProbLinIStateK = matrix.get(x, i, k) + omegaGain - omegaLoss;
                    //if (Double.isNaN(newProbLinIStateK)) {
                        //System.out.println();
                    //}
                    newLineageProbs.add(newProbLinIStateK);
                }
            
                for (int k = 0; k < states; k++) {
                    matrix.set(x, i, k, newLineageProbs.get(k));
                }   
            }
        }
    }
    
    //Set parent lineage state probs after coalescent event
    public void updateProbsAfterCoal(int x, int coalNode, int daughterLineage1Index, int daughterLineage2Index, double lambdaSum, DoubleArrayList lambdaSumForEachK) 
    {
        
        for (int k = 0; k < states; k++) {
            double pThisState = lambdaSumForEachK.get(k) / lambdaSum;
            if (Double.isNaN(pThisState)) {
                //System.out.println("WARNING: Lineage state probs returned NaN after coalescent event");
                pThisState = 0.0;
            }
            if (pThisState < 0) {
                //System.out.println("WARNING: Lineage state probs returned negative after coalescence event");
                pThisState = 0.0;
            }
            matrix.set(x, coalNode, k, pThisState);
        }
    }
    
    public DoubleArrayList getLineagesByState(int x) 
    {
        DoubleArrayList lineagesByState = new DoubleArrayList();
        int numLineages = activeLineages.size();
        for (int k = 0; k < states; k++) {
            double pSum = 0.0;
            for (int lin = 0; lin < numLineages; lin++) {
                int i = activeLineages.get(lin);
                double pMass = matrix.get(x, i, k);
                pSum += pMass;    
            }
            lineagesByState.add(pSum);
        }
        return lineagesByState;
    }
    
}
