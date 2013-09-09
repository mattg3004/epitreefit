package epitreefit;

import java.io.IOException;
import cern.colt.list.*;


/**
 *
 * @author David
 */
public class EpiTreeFit {

    public static void main(String[] args) throws IOException
    {
        
        //Load genealogy
        String treeFileName = "DENV1_subMixed_1130_Fit.tre";
        TreeNode[] tree = NewickReader.getNewickTree(treeFileName);
        
        //Set output file names
        String outputFileStem = "DENV1_subMixed1130"; // + Integer.toString(run+1);
        String outputFilePostFix = "test_082513";
        
        //Reset tip labels to all be in HCMC for unstructured models
        DoubleArrayList tipTimes = TreeUtils.getTipTimes(tree);
        int tips = tipTimes.size();
        for (int tip = 0; tip < tips; tip++) {
            //System.out.println(startTree[tip].nodeLabel);
            tree[tip].nodeLabel = 1;
        }
         
	//Set start and end times and dt integration step
        //double = Double.NEGATIVE_INFINITY; //If NEGATIVE_INFINITY, startTime is root time 
        double startTime = 731582.0; //(Jan_01_2003) for DENV1 tree
	double endTime = Double.POSITIVE_INFINITY; //If POSITIVE_INFINITY, endTime is terminal sampling time
	double dt = 0.5; //Integration dt time step for simulation (in days)
        
        //Specify epidemiological model
        EpiModel epi = new EpiModel2PopDengue(); //for direct transmission models (spatial and non-spatial)
        //EpiModel epi = new EpiModelDengueVectored(); //for vector-borne models (spatial and non-spatial)
        epi.setInitParams();
        
        //Specify coalescent model (needs to match EpiModel)
        StructCoalModel coal = new StructCoalModel2PopDengue(); //for vector-borne models
        //StructCoalModel coal = new StructCoalModelDengueVectored(); //for vector-borne models
        coal.make();
	
	//Set proposal density covariance matrix
        ProposalDensity thetaProposal = new ProposalDensity();
        thetaProposal.setRandomGenerator();
        double [][] thetaCov = {{0.000001,0,0,0},{0,0.0000001,0,0},{0,0,0.0000001,0},{0,0,0,0.000001}}; //for unstructured DENV model w/ four params
        thetaProposal.setCov(thetaCov, 4);
	
        //MCMC params
	int particles = 1; //keep at one (this is a holdover from the particle filtering algorithms in Rasmussen et al. 2011)
	int iterations = 10000; //MCMC iterations
        int outputFreq = 10; //frequency at which MCMC samples are written to file
        boolean verbose = true; //Print out MCMC output?
        
        MetropolisHastings mcmc = new MetropolisHastings();
	mcmc.runMCMC(iterations, particles, outputFreq, tree, startTime, endTime, dt, thetaProposal, epi, coal, outputFileStem, outputFilePostFix, verbose);
        
    }
}
