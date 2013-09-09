package epitreefit;

import cern.colt.list.*;
import cern.colt.matrix.*;
import java.util.ArrayList;
import java.util.Collections;

/**
 * Weights to compute coalescent likelihoods
 * @author David
 */

public class Weights {
    
    DoubleMatrix2D matrix;
    int jParticles;
    int times;
    
    //These are here just to speed things up
    DoubleFactory2D factory2D = DoubleFactory2D.dense;
    //DoubleMatrix2D birthRateMatrix; DoubleMatrix2D migrationRateMatrix; DoubleMatrix2D transMatrix;
    
    public void getMatrix(int jNum, int ts, DoubleFactory2D factory2D)
    {
        jParticles = jNum;
        times = ts- 1;
	matrix = factory2D.make(jParticles, (times));
        
    }
    
    //Compute likelihood over an interval without a coalescent event
    public void updateWeightsNoCoal(DoubleMatrix2D xCurr, LineageStateProbs stateProbs, DoubleArrayList theta, int time, StructCoalModel coal, double dtTime, double currTime)
    {

        int states = coal.Y.size();
        int numLineages = stateProbs.activeLineages.size();
        double lambdaSum; int i; int j; double popCoalRate; double pIIsInK; double pIIsInL; double pJIsInK; double pJIsInL; double linPairCoalRate;
        double Ak; double Al; //double rateAkAk; double rateAkAl; double rateAlAl; double lambdaSumAg;
        DoubleArrayList linesAk = new DoubleArrayList();
        
        if (numLineages > 1) {

            for (int x = 0; x < jParticles; x++) { //use x's for particle indexes to avoid confusion with lineage indexes i and j

                double pWeight = 1.0;
                
                //Update F matrix and G matrix and population sizes in Y
                coal.updateF(xCurr.viewColumn(x), theta, currTime);
                coal.updateG(xCurr.viewColumn(x), theta);
                coal.updateY(xCurr.viewColumn(x));
                
                linesAk = stateProbs.getLineagesByState(x); //returns number of lineages in each state k 
                
                boolean sizeFlag = false;
                boolean lineFlag = false;
                for (int st = 0; st < states; st++) {
                    if (linesAk.get(st) > coal.Y.getQuick(st) || Double.isNaN(linesAk.get(st))) {
                        sizeFlag = true; //Population size is incompatible with number of lineages in state 
                        
                    }
                    if (linesAk.get(st) < 1.0) {
                        lineFlag = false; //Fewer than one lineage in state
                    }
                }
                
                if (sizeFlag) {
                    pWeight = 0.0; //set likelihood equal to zero
                } else {
                
                    lambdaSum = 0.0; //total rate of coalescence across all pairs of lineages
                    
                    //Fast method if all A_k greater than 1.0 for all states k (groups lineages by states)
                    if (lineFlag == false)  {
                        
                        //General case for m > 1
                        for (int k = 0; k < states; k++) {
                            if (coal.Y.getQuick(k) > 0.0) {
                                Ak = linesAk.getQuick(k);
                                lambdaSum += ((Ak*(Ak-1))/2) * (2 * coal.F.getQuick(k,k) / (coal.Y.getQuick(k)*coal.Y.getQuick(k)));
                            }
                        }
                        
                        for (int k = 0; k < states; k++) {
                            for (int l = 0; l < k; l++) { //l < k so not double counting here
                                if (coal.Y.getQuick(k) > 0.0 & coal.Y.getQuick(l) > 0.0) {
                                    Ak = linesAk.getQuick(k);
                                    Al = linesAk.getQuick(l);
                                    lambdaSum += (Ak * Al) * ((coal.F.getQuick(k,l) + coal.F.getQuick(l,k)) / (coal.Y.getQuick(k)*coal.Y.getQuick(l)));
                                }
                            }
                        }
                     
                    //Slow method for computing lambda for each pair of lineages (see Volz, Genetics, 2012)    
                    } else {

                        for (int linI = 0; linI < numLineages; linI++) { //for all lineages i

                            i = stateProbs.activeLineages.get(linI);

                            for (int linJ = linI+1; linJ < numLineages; linJ++) { //for all unique pairs of lineages i and j

                                j = stateProbs.activeLineages.get(linJ);

                                //Compute probability of lineages i and j NOT coalescing
                                for (int k = 0; k < states; k++) { //for all states k
                                    if (coal.Y.getQuick(k) <= 0.0) {
                                        //lambdaSum can not increase
                                    } else {
                                        for (int l = 0; l < states; l++) {//for all state l
                                            if (coal.Y.getQuick(l) <= 0.0) {
                                                //lambdaSum can not increase
                                            } else {
                                                //popCoalRate = coal.computeCoalRateOneDirection(k, l);
                                                popCoalRate = coal.F.getQuick(k,l) / (coal.Y.getQuick(k)*coal.Y.getQuick(l));
                                                pIIsInK = stateProbs.matrix.getQuick(x, i, k); //(particle, lineage, state)
                                                pIIsInL = stateProbs.matrix.getQuick(x, i, l);
                                                pJIsInK = stateProbs.matrix.getQuick(x, j, k);
                                                pJIsInL = stateProbs.matrix.getQuick(x, j, l);
                                                linPairCoalRate = popCoalRate * (pIIsInK*pJIsInL + pIIsInL*pJIsInK);
                                                if (Double.isNaN(linPairCoalRate)) { //Checking for NaN's here might not be needed
                                                    linPairCoalRate = 0.0;
                                                }
                                                lambdaSum += linPairCoalRate;

                                            }
                                        }
                                    }
                                }
                            } 
                        }
                    }

                    pWeight = Math.exp(-lambdaSum * dtTime); //coalescent likelihood over interval
                }
                if (Double.isNaN(pWeight)) {
                    System.out.println("WARNING: Some particle weights returned NaN");
                    pWeight = Double.MIN_VALUE;
                }
                if (pWeight <= 0.0) {
                    pWeight = Double.MIN_VALUE;
                    //System.out.println("WARNING: Zero or negative particle weights found");
                }
                matrix.setQuick(x,time-1,pWeight);
                
            }
        } else { //less than two lineages
            for (int x = 0; x < jParticles; x++) {
                matrix.setQuick(x,time-1,1.0); //assign weight = 1.0 if less than two lineages present
            }
        }

    }
    
    //Compute likeliood of coalescent event
    public void updateWeightsCoal(ZVectors dataZ, int zLoc, TreeNode[] tree, DoubleMatrix2D xCurr, LineageStateProbs stateProbs, DoubleArrayList theta, int time, StructCoalModel coal, double currTime) 
    {
    
        int totalEvents = dataZ.omegaEvents.get(zLoc).size();
        int states = coal.Y.size();
        
        //Establish precedence among coalescing lineages
        ArrayList<Integer> coalNodesAtEvent = new ArrayList<Integer>();
        for (int event = 0; event < totalEvents; event++) {
            if (dataZ.omegaEvents.get(zLoc).get(event)==1) {
                int coalNode = dataZ.nodePointers.get(zLoc).get(event);
                coalNodesAtEvent.add(coalNode);
            }
        }
        Collections.sort(coalNodesAtEvent);
        Collections.reverse(coalNodesAtEvent);
        
        int coalEvents = coalNodesAtEvent.size();
        for (int event = 0; event < coalEvents; event++) {
                
            int coalNode = coalNodesAtEvent.get(event);
            int daughterLineage1Index = tree[coalNode].childNodes[0] - 1;
            int daughterLineage2Index = tree[coalNode].childNodes[1] - 1;

            for (int x = 0; x < jParticles; x++) { //use x's for particle indexes to avoid confusion with lineage indexes i and j

                coal.updateF(xCurr.viewColumn(x), theta, currTime);
                coal.updateG(xCurr.viewColumn(x), theta);
                coal.updateY(xCurr.viewColumn(x));

                double lambdaSum = 0.0;
                DoubleArrayList lambdaSumForEachK = new DoubleArrayList();
                for (int k = 0; k < states; k++) { //for all states k
                    double lambdaSumK = 0.0;
                    if (coal.Y.getQuick(k) <= 0.0) {
                        //lambdaSumK does not increase
                    } else {
                        for (int l = 0; l < states; l++) {//for all state l
                            if (coal.Y.getQuick(l) <= 0.0) {
                                //lambdaSum can not increase
                            } else {
                                double popCoalRate = coal.computeCoalRateOneDirection(k, l);
                                double pIIsInK = stateProbs.matrix.getQuick(x, daughterLineage1Index, k); //(particle, lineage, state)
                                double pIIsInL = stateProbs.matrix.getQuick(x, daughterLineage1Index, l);
                                double pJIsInK = stateProbs.matrix.getQuick(x, daughterLineage2Index, k);
                                double pJIsInL = stateProbs.matrix.getQuick(x, daughterLineage2Index, l);
                                double linPairCoalRate = popCoalRate * (pIIsInK*pJIsInL + pIIsInL*pJIsInK);
                                if (Double.isNaN(linPairCoalRate)) { //checking for NaN's here might not be needed
                                    linPairCoalRate = 0.0;
                                }
                                if (Double.isInfinite(linPairCoalRate)) {
                                    linPairCoalRate = 0.0;
                                }
                                lambdaSum += linPairCoalRate;
                                lambdaSumK += linPairCoalRate;
                            }
                        }
                    }
                    lambdaSumForEachK.add(lambdaSumK);
                }
                
                //if (lambdaSum == 0) {
                    //System.out.println("WARNING: Found zero probability of a coalescence event");
                //}

                //Update weight to reflect coalescent event}
                double pWeightNew = matrix.get(x,time-1) * lambdaSum;
                if (Double.isNaN(pWeightNew)) {
                    System.out.println("WARNING: Some particle weights returned NaN");
                }
                if (pWeightNew <= 0.0) {
                    pWeightNew = Double.MIN_VALUE;
                    //System.out.println("WARNING: Zero  or negative particle weights found");
                }
                matrix.setQuick(x,time-1,pWeightNew);

                //Update LineageStateProbs
                stateProbs.updateProbsAfterCoal(x, coalNode, daughterLineage1Index, daughterLineage2Index, lambdaSum, lambdaSumForEachK);
            }

            //Remove daughter lineages from and add parent lineage to activeLineages
            int listIndex1 = stateProbs.activeLineages.indexOf(daughterLineage1Index);
            if (listIndex1 != -1) {
                stateProbs.activeLineages.remove(listIndex1);
            } else {
                System.out.println("Cannot remove daughter lineage");
            }

            int listIndex2 = stateProbs.activeLineages.indexOf(daughterLineage2Index);
            if (listIndex2 != -1) {
                stateProbs.activeLineages.remove(listIndex2);
            } else {
                System.out.println("Cannot remove daughter lineage");
            }

            stateProbs.lineagesRemoved.add(daughterLineage1Index);
            stateProbs.lineagesRemoved.add(daughterLineage2Index);

            stateProbs.activeLineages.add(coalNode);
            stateProbs.lineagesAdded.add(coalNode);     
            
        }
    }
    
    public double computeLogLikelihood()
    {
	double totalLogLike = 0;
        double sumOfWeights;
	double avgWeight;
	for (int n = 0; n < times; ++n) {
            sumOfWeights = 0;
            for (int j = 0; j < jParticles; ++j) { //only summing over one particle here for deterministic simulations
                sumOfWeights += matrix.getQuick(j,n);
            }
            avgWeight = sumOfWeights / jParticles;
            totalLogLike += Math.log(avgWeight);
	}	
	return totalLogLike;
        
    } //END method
    
}

