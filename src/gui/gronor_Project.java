package gronor;

import java.util.*;
import java.util.StringTokenizer;
import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.table.*;
import javax.swing.event.*;
import javax.swing.border.*;
import ptolemy.plot.*;

import cformat.PrintfWriter;

public class gronor_Project extends JFrame implements ActionListener, ChangeListener, WindowListener, MouseListener, ItemListener{
	
	private static final long serialVersionUID = 4L;

	JMenuBar menubar = new JMenuBar();
	JMenuItem create, open, close, list, save, quit;
	String projectName, projectFile;

	JPanel parametersPanel;
	JPanel optionsPanel;
	JPanel dimensionPanel;
	JPanel jobPanel;
	JPanel expansionPanel;
	JPanel numberPanel;
	JPanel threshPanel;
	JPanel basisPanel;
	JPanel statesPanel;
	JPanel fragmentsPanel;
	JPanel fragmentsButtonsPanel;
	JPanel energiesPanel;
	JPanel mebfsPanel;
	JPanel nociResultsPanel;
	JPanel nociEnergiesPanel;
	JPanel nociPlotPanel;

	Integer maxSets = 6;	// Maximum number of state definitions
	
	Integer maxFragments = 32;	// Maximum number of fragments
	
	Integer numStateNames = 15;
	
	Integer numSets=0;			// Number of state sets
	
	Integer numFragments;		// Number of fragments
	Integer numFragmentStates=0;  // Number of unique fragment states used

	Integer[] lenStateList = new Integer[maxSets];				// Number of states included from list for state set
	Integer[] spinStateList = new Integer[maxSets];				// Spin of states included from list for state set
	Integer[][] ndxStateList = new Integer[maxSets][15];		// Index to stateNames[] for each state in state set
	Integer[] fragmentStates = new Integer[35];					// Index to stateNames[] for each unique fragment state used
	
	Integer[] ndxMebfTable = new Integer[35];					// Index of each state into into MEBF table 
	Integer[] ndxMebfState = new Integer[35];					// Index into MEBF table for each state in MEBF 
	
	Integer[][] dimFragments = new Integer[maxFragments][11];	// Dimensions of fragments: number of atoms, states, electrons, etc.
	String[] namFragments = new String[maxFragments];			// Names of fragments used to determine xyz-formatted coordinate origin 
	Double[][] movFragments = new Double[maxFragments][6];		// Rotation and translation of coordinates with respect to original source	
	Double[] energiesDFT = new Double[maxFragments];			// DFT optimized energy of S0 state from NWChem
	Double[] energiesSCF = new Double[maxFragments];			// SCF energy of S0 state from Molcas
	Double[][] energiesCASSCF = new Double[maxFragments][12];	// CASSCF energies of all states of fragment from Molcas
	Double[][] energiesCASPT2 = new Double[maxFragments][12];	// CASPT2 energies of all states of fragment from Molcas
	Object[][] dimensionData;
	
	Integer numRanks=12;		// Number of ranks in internal mpirun
	Integer numThreads=12;
	Integer memory=2048;
	Integer expansion=0;
	Integer basisSet=0;
	Integer contract=1;
	Integer cholesky=2;
	
	String account="chm154";
	String jobName="GronOR";
	String timeLimit="00:30:00";
	
	Integer numEnergies = 0;	// Number of energy entries

	Integer numMEBFRows = 0;
    Integer numMEBFs = 0;
    Integer newMEBFs = 0;
    Integer maxMEBFs = 36;
    Integer maxMer = 5;
    Integer maxMEBFstates=34;
    
    Integer numStateLabels = 12;

	Integer[][] numEnergiesGronOR = new Integer[maxMEBFs][2];		// Number of GronORenergy entries
	Double[][][] energiesGronOR = new Double[35][maxMEBFs][2];
	Double[][][] energiesRelGronOR = new Double[35][maxMEBFs][2];
	Double[][][] coefGronOR = new Double[35][35][maxMEBFs];

    String[] stateListLabels = new String[] {"ID","Num","Spin","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"};
    String[] fragmentLabels = new String[] {"ID", "SRC", "XYZ File", "Atoms", "States", "Electrons", "CASe", "CASo", "Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "Alt"};
    String[] fragmentNames = new String[] {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"};
    String[] energyNames = new String[] {"E(DFT)", "E(SCF)", "E(CASSCF)", "E(CASPT2)"};
    String[] stateLabels = new String[] {" ", "S0", "S1", "T1", "D-", "D+", "S2", "T2", "E-", "E+"," "," "};

    String[] stateNames = new String[] {"S0","S1","S2","D0","D1","T1","T2","S+","D+","T+","S-","D-","T-","Q1"};
    Integer[] stateSpins = new Integer[] {1,1,1,2,2,3,3,1,2,3,1,2,3,5};
	String[] mebfLabels = new String[] {"ID","n-mer","Spin","Charge","States","Frag","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36"};
	String[] nociLabels = new String[] {"ID", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9" };

	String[] expMEBF = new String[] {"small", "medium", "high"};
	String[] basisSets = new String[] {"ano-s", "ano-l", "ano-ccr", "ano-ccr-vdzp"};
	String[] contracts = new String[] {"s", "m","l"};
	String[] choleskys = new String[] {"high", "medium", "low", "none"};
	
	Integer[] stateLabelIndex = new Integer[numStateLabels+3];
	
    JTable statesTable;
    JTable fragmentsTable;
    JTable energiesTable;
    JTable mebfsTable;
	JTable dimensionTable = new JTable();
	JTable jobTable = new JTable();
	JTable numberTable = new JTable();
	JTable threshTable = new JTable();
	JTable basisTable = new JTable();
	JTable expansionTable = new JTable();
	JTable contractTable = new JTable();
	JTable choleskyTable = new JTable();
    JTable nociEnergiesTable;
	DefaultTableModel nociEnergiesTableModel = new DefaultTableModel();
    
    Boolean numFragmentsChanged = false;

	DefaultTableModel statesTableModel = new DefaultTableModel(null,stateListLabels);
	DefaultTableModel fragmentsTableModel = new DefaultTableModel(null,fragmentLabels);
	DefaultTableModel energiesTableModel = new DefaultTableModel(null,stateLabels);
	DefaultTableModel mebfsTableModel = new DefaultTableModel(null,mebfLabels);
	DefaultTableModel mebfsEnergiesTableModel = new DefaultTableModel(null,nociLabels);

    Object[][] stateDefinitions = new Object[maxSets][18];
    Object[][] fragmentDefinitions = new Object[maxFragments][15];
    Object[][] stateEnergies = new Object[4*maxFragments][12];
    Object[][] mebfDefinitions = new Object[maxMEBFs*maxMer][maxMEBFstates+6];
    Object[][] mebfEnergies = new Object[maxMEBFs][10];
    Object[][] nociEnergies = new Object[maxMEBFstates][maxMEBFs];
    
    Integer[][] mebfIndex = new Integer[maxMEBFs*maxMer][2];
    Integer[][] mebfSpecification = new Integer[maxMEBFs][5];	 	        	// 0: n-mer; 1: spin; 2: charge; 3:number of states to include 4:fragment; 
    String[] mebfName = new String[maxMEBFs];							        // name of mefb, e.g. AB for fragment A,B combination
    Integer[][][] mebfFragments = new Integer[maxMEBFs][maxMer][maxMEBFstates];	// index to fragment states included in this mebf
    
    Integer numberStateEnergies = 0;
    Integer numberMebfDefinitions = 0;
    Integer numberMebfEnergies = 0;
    Integer numberFragmentDefinitions = 0;
    Integer numberStateDefinitions = 0;
    
    Integer[][] energyFragment = new Integer[maxFragments][2];

    gronor_Fragment fragment = new gronor_Fragment();

    Double[] RandT = new Double[6];
    
    Double thresh_CI = 1.0e-5;
    Double thresh_MO = 1.0e-4;
    Double ipea = 0.25;

    Integer noci0=0;
    Integer noci1=0;

    Integer[][][][] couple3 = new Integer[5][5][5][5];		        // indices: target spin; frag 1,2,3 spin;
    Integer[][][][][][] couple4 = new Integer[5][5][5][5][5][2];	// indices: target spin; frag 1,2,3,4 spin; coupling 1,2
    
    JFrame dialogFrame;

    JFileChooser chooser;
    extensionFilter filter;
    String prjFile;
        
    Integer numAtoms=12, numStates=0, newFragments=0, newSets=0;
   
    String[] methods = new String[] {"SCF", "CASSCF", "CASPT2"};
    JLabel[] methodLabel = new JLabel[3];
    
    Integer[] fragments;
	
    Integer numDFT=0, numSCF=0, numCASSCF=0, numCASPT2=0;
    
    Integer numCASe=8, numCASo=8, numAlt=0, numOcc=0;
    Integer[][] alter = new Integer[12][2];
    
    Plot energyDiagram = new Plot();
    
    Container container = null;
    
	gronor_Project(String pName){
		
		super("GronOR Project "+pName);
	    super.setSize(1200,900);
	    super.setLocation(150,150);
	    super.setBackground(Color.lightGray);
	    super.setForeground(Color.blue);
	    super.addWindowListener(this);
	    super.setFont(new Font("Serif",Font.BOLD,10));

		projectName=pName;
		projectFile=projectName+".prj";
		
	    this.setJMenuBar(menubar);

	    JMenu select = new JMenu("File");
	    menubar.add(select);
	    
	    select.add(save = new JMenuItem("Save"));
	    save.addActionListener(new ActionListener(){
	      public void actionPerformed(ActionEvent e){ 
	    	  if(!writeProjectFile(projectFile)) System.exit(0);
	    	  writeClearScripts();
	    }});


		if(!readProjectFile(projectFile)) initializeProject();
		
		initializeWindow();
		
	    setVisible(true);
	    
	}

	private void writeClearScripts() {
		String fileName = "clear.run";
		try {
			PrintfWriter clearFile = new PrintfWriter(new FileWriter(fileName));
			clearFile.println("#!/usr/bin/tcsh");
			clearFile.println("rm -f *.LprOrb");
			clearFile.println("rm -f *.oneint");
			clearFile.println("rm -f *.RasOrb*");
			clearFile.println("rm -f INPORB");
			clearFile.println("rm -f *.INPORB*");
			clearFile.println("rm -f *.ROTORB*");
			clearFile.println("rm -f *.output");
			clearFile.println("rm -f *.RUNFIL");
			clearFile.println("rm -f *.RUNFILE");
			clearFile.println("rm -f *.VecDet*");
			clearFile.println("rm -f *.molden*");
			clearFile.println("rm -f *.det");
			clearFile.println("rm -f *.vec");
			clearFile.println("rm -f *.status");
			clearFile.println("rm -f *.one");
			clearFile.println("rm -f *.two");
			clearFile.println("rm -f *.sys");
			clearFile.println("rm -f *.nwout");
			clearFile.println("rm -f *.hess");
			clearFile.println("rm -f *.TRAONE");
			clearFile.println("rm -f *.c");
			clearFile.println("rm -f *.db");
			clearFile.println("rm -f *.movecs");
			clearFile.println("rm -f *.p");
			clearFile.println("rm -f CH*");
			clearFile.println("rm -f _CHMOT*");
			clearFile.println("rm -f COMMON*");
			clearFile.println("rm -f *.GssOrb");
			clearFile.println("rm -f *_CHMOT*");
			clearFile.println("rm -f *.CH*");
			clearFile.println("rm -f *.COMMON*");
			clearFile.println("rm -f *.ScfOrb");
			clearFile.println("rm -f *.SpdOrb*");
			clearFile.println("rm -f *.zmat");
			clearFile.println("rm -f *.b*");
			clearFile.println("rm -f ONEINT*");
			clearFile.println("rm -f RUNFIL*");
			clearFile.println("rm -f TRA*");
			clearFile.println("rm -f xml*");
			clearFile.println("rm -f *.arx");
			clearFile.println("rm -f *.lus");
			clearFile.println("rm -f *.out");
			clearFile.println("rm -f *.day");
			clearFile.println("rm -f *.cml");
			clearFile.println("rm -f *.pro");
			clearFile.println("rm -f *.rnk");
			clearFile.println("rm -f *.tst");
			clearFile.println("rm -f *.cpr");
			clearFile.println("rm -f *.xmldump");
			for(int i=0; i<numFragments; i++) {
				clearFile.println("rm -f "+projectName.trim()+fragmentDefinitions[i][0]+".xyz");
			}
			for(int i=0; i<numMEBFs; i++) {
				clearFile.println("rm -f "+projectName.trim()+mebfName[i]+".xyz");
			}
			clearFile.close();
		} catch(IOException ei) {
		}
		fileName = "clearall.run";
		try {
			PrintfWriter clearFile = new PrintfWriter(new FileWriter(fileName));
			clearFile.println("#!/usr/bin/tcsh");
			clearFile.println("rm -f "+projectName.trim()+"*.input");
			clearFile.println("rm -f "+projectName.trim()+"*.inp");
			clearFile.println("rm -f "+projectName.trim()+"*.nw");
			clearFile.println("rm -f "+projectName.trim()+"*.run");
			clearFile.println("rm -f "+projectName.trim()+"*.lsf");
			clearFile.println("rm -f "+projectName.trim()+"*.slurm");
			clearFile.println("rm -f *.LprOrb");
			clearFile.println("rm -f *.oneint");
			clearFile.println("rm -f *.RasOrb*");
			clearFile.println("rm -f INPORB");
			clearFile.println("rm -f *.INPORB*");
			clearFile.println("rm -f *.ROTORB*");
			clearFile.println("rm -f *.output");
			clearFile.println("rm -f *.RUNFIL");
			clearFile.println("rm -f *.RUNFILE");
			clearFile.println("rm -f *.VecDet*");
			clearFile.println("rm -f *.molden*");
			clearFile.println("rm -f *.det");
			clearFile.println("rm -f *.vec");
			clearFile.println("rm -f *.status");
			clearFile.println("rm -f *.one");
			clearFile.println("rm -f *.two");
			clearFile.println("rm -f *.sys");
			clearFile.println("rm -f *.nwout");
			clearFile.println("rm -f *.hess");
			clearFile.println("rm -f *.TRAONE");
			clearFile.println("rm -f *.c");
			clearFile.println("rm -f *.db");
			clearFile.println("rm -f *.movecs");
			clearFile.println("rm -f *.p");
			clearFile.println("rm -f CH*");
			clearFile.println("rm -f _CHMOT*");
			clearFile.println("rm -f COMMON*");
			clearFile.println("rm -f *.GssOrb");
			clearFile.println("rm -f *_CHMOT*");
			clearFile.println("rm -f *.CH*");
			clearFile.println("rm -f *.COMMON*");
			clearFile.println("rm -f *.ScfOrb");
			clearFile.println("rm -f *.SpdOrb*");
			clearFile.println("rm -f *.zmat");
			clearFile.println("rm -f *.b*");
			clearFile.println("rm -f ONEINT*");
			clearFile.println("rm -f RUNFIL*");
			clearFile.println("rm -f TRA*");
			clearFile.println("rm -f xml*");
			clearFile.println("rm -f *.arx");
			clearFile.println("rm -f *.lus");
			clearFile.println("rm -f *.out");
			clearFile.println("rm -f *.day");
			clearFile.println("rm -f *.cml");
			clearFile.println("rm -f *.pro");
			clearFile.println("rm -f *.rnk");
			clearFile.println("rm -f *.tst");
			clearFile.println("rm -f *.cpr");
			clearFile.println("rm -f *.xmldump");
			clearFile.println("rm -f *.orb");
			clearFile.println("rm -f *.runfil");
			clearFile.println("rm -f *.comorb");
			clearFile.println("rm -f *.ONEINT");
			for(int i=0; i<numFragments; i++) {
				clearFile.println("rm -f "+projectName.trim()+fragmentDefinitions[i][0]+".xyz");
			}
			for(int i=0; i<numMEBFs; i++) {
				clearFile.println("rm -f "+projectName.trim()+mebfName[i]+".xyz");
			}
			clearFile.close();
		} catch(IOException ei) {
		}
	}
	
	private Boolean readProjectFile(String fileName) {
		try {
		    	BufferedReader br = new BufferedReader(new FileReader(fileName));
		    	String card;
			    card=br.readLine();
			    jobName=card.trim();
			    card=br.readLine();
			    account=card.trim();
			    card=br.readLine();
			    timeLimit=card.trim();
		    	card=br.readLine();
			    numSets = Integer.valueOf(card.substring(0,6).trim());
		    	newSets = numSets;
			    numFragments = Integer.valueOf(card.substring(6,12).trim());
		    	newFragments = numFragments;
		    	numMEBFs = Integer.valueOf(card.substring(12,18).trim());
		    	newMEBFs = numMEBFs;
		    	numRanks = Integer.valueOf(card.substring(18,24).trim());
		    	numThreads = Integer.valueOf(card.substring(24,30).trim());
		    	basisSet  = Integer.valueOf(card.substring(30,36).trim());
		    	contract  = Integer.valueOf(card.substring(36,42).trim());
		    	cholesky  = Integer.valueOf(card.substring(42,48).trim());
		    	expansion = Integer.valueOf(card.substring(48,54).trim());
		    	thresh_MO = Double.valueOf(card.substring(54,74)).doubleValue();
		    	thresh_CI = Double.valueOf(card.substring(74,94)).doubleValue();
		    	ipea      = Double.valueOf(card.substring(94,100)).doubleValue();
			    
		    	for(int i=0; i<numSets; i++) {
		    		card=br.readLine();
		    		lenStateList[i]=Integer.valueOf(card.substring(0,6).trim());
		    		spinStateList[i]=Integer.valueOf(card.substring(6,12).trim());
		    		for(int j=0; j<lenStateList[i]; j++) {
		    			ndxStateList[i][j]=Integer.valueOf(card.substring((j+2)*6,(j+3)*6).trim());
		    		}
		    	}
			    for(int i=0; i<numFragments; i++) {
			    	card=br.readLine();
			    	namFragments[i]=card.trim();
			    	card=br.readLine();
			    	for(int j=0; j<6; j++) movFragments[i][j]=Double.valueOf(card.substring(j*20,j*20+20)).doubleValue();
			    	card=br.readLine();
			    	dimFragments[i][0]=Integer.valueOf(card.substring(0,6).trim());    // Index to XYZ
			    	dimFragments[i][1]=Integer.valueOf(card.substring(6,12).trim());   // Number of atoms
			    	dimFragments[i][2]=Integer.valueOf(card.substring(12,18).trim());  // Index into states list
			    	dimFragments[i][3]=Integer.valueOf(card.substring(18,24).trim());  // Number of electrons
			    	dimFragments[i][4]=Integer.valueOf(card.substring(24,30).trim());  // Number of CAS electrons
			    	dimFragments[i][5]=Integer.valueOf(card.substring(30,36).trim());  // Number of CAS orbitals
			    	dimFragments[i][6]=Integer.valueOf(card.substring(36,42).trim());  // Number of DFT energy
			    	dimFragments[i][7]=Integer.valueOf(card.substring(42,48).trim());  // Number of SCF energy
			    	dimFragments[i][8]=Integer.valueOf(card.substring(48,54).trim());  // Number of CASSCF energies
			    	dimFragments[i][9]=Integer.valueOf(card.substring(54,60).trim());  // Number of CASPT2 energies
			    	dimFragments[i][10]=Integer.valueOf(card.substring(60,66).trim()); // Number of orbital alterations
			    	fragmentDefinitions[i][1]=card.substring(71,72).trim(); // Source fragment number
			    	card=br.readLine();
//			    	if(dimFragments[i][6]>0) energiesDFT[i]=Double.valueOf(card.substring(0,20)).doubleValue();
//			    	if(dimFragments[i][7]>0) energiesSCF[i]=Double.valueOf(card.substring(20,40)).doubleValue();
			    	energiesDFT[i]=Double.valueOf(card.substring(0,20)).doubleValue();
			    	energiesSCF[i]=Double.valueOf(card.substring(20,40)).doubleValue();
			    	card=br.readLine();
			    	for(int j=0; j<dimFragments[i][8]; j++) energiesCASSCF[i][j]=Double.valueOf(card.substring(j*20,j*20+20)).doubleValue();
			    	card=br.readLine();
			    	for(int j=0; j<dimFragments[i][9]; j++) energiesCASPT2[i][j]=Double.valueOf(card.substring(j*20,j*20+20)).doubleValue();
			    }
			    
				for(int i=0; i<numMEBFs; i++) {
			    	card=br.readLine();
					mebfName[i]=card.trim();
			    	card=br.readLine();
					for(int j=0; j<5; j++) mebfSpecification[i][j]=Integer.valueOf(card.substring(j*6,j*6+6).trim());
					Integer nmer = mebfSpecification[i][0];
					for(int j=0; j<nmer; j++) {
				    	card=br.readLine();
						for(int k=0; k<maxMEBFstates; k++) mebfFragments[i][j][k]=Integer.valueOf(card.substring(k*6,k*6+6).trim());
					}
				}
				
				br.close();
				return true;
			} catch(Exception ee) {
				return false;
			}
	}

	private Boolean writeProjectFile(String fileName) {
		numFragmentsChanged = numFragments!=newFragments;
		try {
			    PrintfWriter fw = new PrintfWriter(new FileWriter(fileName));
			    if(newFragments<numFragments) numFragments=newFragments;
			    fw.println(jobName.trim());
			    fw.println(account.trim());
			    fw.println(timeLimit.trim());
			    fw.printf("%6d",newSets);
			    fw.printf("%6d",newFragments);
			    fw.printf("%6d",numMEBFs);
			    fw.printf("%6d",numRanks);
			    fw.printf("%6d",numThreads);
			    fw.printf("%6d",basisSet);
			    fw.printf("%6d",contract);
			    fw.printf("%6d",cholesky);
			    fw.printf("%6d",expansion);
			    fw.printf("%20.10f",thresh_MO);
			    fw.printf("%20.10f",thresh_CI);
			    fw.printf("%6.3f",ipea);
			    fw.println();
		    	for(int i=0; i<numSets; i++) {
		    		fw.printf("%6d",lenStateList[i]);
		    		fw.printf("%6d",spinStateList[i]);
		    		for(int j=0; j<lenStateList[i]; j++) {
			    		fw.printf("%6d",ndxStateList[i][j]);
		    		}
		    		fw.println();
		    	}
				for(int i=0; i<numFragments; i++) {
					fw.println(namFragments[i].trim());
					for(int j=0; j<6; j++) fw.printf("%20.10f",movFragments[i][j]);
				    fw.println();
				    for(int j=0; j<11; j++) fw.printf("%6d",dimFragments[i][j]);
				    String name = (String) fragmentDefinitions[i][1];
				    fw.printf("%s","     "); fw.printf("%s",name.trim());
				    fw.println();
				    if(dimFragments[i][6]>0) {
				    	fw.printf("%20.10f",energiesDFT[i]);
				    } else {
				    	fw.printf("%20.10f",0.0);
				    }
				    if(dimFragments[i][7]>0) {
				    	fw.printf("%20.10f",energiesSCF[i]);
				    }else {
				    	fw.printf("%20.10f",0.0);
				    }
				    fw.println();
				    for(int j=0; j<dimFragments[i][8]; j++) fw.printf("%20.10f",energiesCASSCF[i][j]);
				    fw.println();
				    for(int j=0; j<dimFragments[i][9]; j++) fw.printf("%20.10f",energiesCASPT2[i][j]);
				    fw.println();
				}
				if(newFragments>numFragments) {
					for(int i=numFragments; i<newFragments; i++) {
						fw.println();
					for(int j=0; j<6; j++) fw.printf("%20.10f",0.0);
				    fw.println();
				    for(int j=0; j<2; j++) fw.printf("%6d",0); 
				    fw.printf("%6d",5); 
				    fw.printf("%6d",0); 
				    fw.printf("%6d",8); 
				    fw.printf("%6d",8);
				    for(int j=6; j<10; j++) fw.printf("%6d",0);
				    fw.println();
				    fw.println();
				    fw.println();
				    fw.println();
					}
				}
				
				for(int i=0; i<numMEBFs; i++) {
					fw.println(mebfName[i].trim());
					for(int j=0; j<5; j++) fw.printf("%6d",mebfSpecification[i][j]); fw.println();
					Integer nmer = mebfSpecification[i][0];
					for(int j=0; j<nmer; j++) {
						for(int k=0; k<maxMEBFstates; k++) fw.printf("%6d",mebfFragments[i][j][k]); fw.println();
					}
				}
				
				fw.close();
				numFragments=newFragments;
				return true;
			} catch(Exception ee) {
				System.out.println("Project file could not be written");
				return false;
			}
	}
	
	private void initializeProject() {
		
		numSets=1;
		newSets=1;
		numFragments=0;
		newFragments=0;
	    numMEBFs = 0;
	    newMEBFs = 0;

		// Initialize with two state sets for S0 and D0 ground states
	    
		for(int i=0; i<maxSets; i++) {
			lenStateList[i]=0;
			spinStateList[i]=1;
			for(int j=0; j<15; j++) ndxStateList[i][j]=0;

			lenStateList[0]=5;
			ndxStateList[0][0]=0;
			ndxStateList[0][1]=1;
			ndxStateList[0][2]=5;
			ndxStateList[0][3]=8;
			ndxStateList[0][4]=11;
			
//			lenStateList[1]=6;
//			ndxStateList[1][0]=3;
//			ndxStateList[1][1]=4;
//			ndxStateList[1][2]=7;
//			ndxStateList[1][3]=9;
//			ndxStateList[1][4]=10;
//			ndxStateList[1][5]=12;
		}
	}
	
	private void createFragments() {
		
		fragmentDefinitions = new Object[numFragments][15];
		
		for(int i=0; i<numFragments; i++) fragmentDefinitions[i][0]=fragmentNames[i];
		
	}
	
	private String getFileRoot(String extension) {
		JFileChooser chooser;
	    extensionFilter filter;
	    String extFile;
	    String rootName;
	    JFrame dialogFrame;
		File prjF = new File("./");
	    FilenameFilter fnFilter = new FilenameFilter() {
	    	public boolean accept(File f, String name) {
	    		return name.endsWith(extension);
	    	}
	    };
	    File[] files = prjF.listFiles(fnFilter);
	    if(files.length==1) {
	    	extFile=files[0].getName();
	    } else {
	    	dialogFrame = new JFrame();
	    	dialogFrame.setSize(300,400);
	    	filter = new extensionFilter(extension);
	    	chooser = new JFileChooser("./");
	    	chooser.setFileFilter(filter);
	    	chooser.showOpenDialog(dialogFrame);
	    	try {
	    	extFile = chooser.getSelectedFile().toString();
	    	} catch (NullPointerException e1) {
	    		return "";
	    	}
	    };
		rootName=extFile.substring(extFile.lastIndexOf("/")+1,extFile.indexOf(extension));
		return rootName;
	}

	public Integer getNumAxyz(String root) {
		String fileName=root.trim()+".xyz";
		String card;
		StringTokenizer st;
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			card=br.readLine();st = new StringTokenizer(card," ");
			Integer nA=Integer.valueOf(st.nextToken());
			br.close();
			return nA;
		} catch(IOException e) {
			System.out.println("IOException in XYZ file "+fileName);
			return 0;
		}
	}

	public Integer getNumHxyz(String root) {
		String fileName=root+".xyz";
		String card;
		StringTokenizer st;
		Integer number = 0;
		Integer nA = 0;
		String aL;
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			card=br.readLine();
			st = new StringTokenizer(card," ");
			nA=Integer.valueOf(st.nextToken());
			card=br.readLine();
			number=0;
			for(int i=0; i<nA; i++) {
				card=br.readLine();
				st = new StringTokenizer(card," ");
				aL=st.nextToken();
				if(aL.equals("H")) number=number+1;
			}
			br.close();
			return number;
		} catch(IOException ef) {
			System.out.println("IOException in XYZ file "+fileName);
			return 0;
		}
	}
	
	public Integer getNumExyz(String root) {
		String fileName=root+".xyz";
		String card;
		StringTokenizer st;
		Integer number = 0;
		Integer nA = 0;
		String aL;
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			card=br.readLine();
			st = new StringTokenizer(card," ");
			nA=Integer.valueOf(st.nextToken());
			card=br.readLine();
			number=0;
			for(int i=0; i<nA; i++) {
				card=br.readLine();
				st = new StringTokenizer(card," ");
				aL=st.nextToken();
				if(aL.equals("H"))  number=number+1;
				if(aL.equals("He")) number=number+2;
				if(aL.equals("Li")) number=number+3;
				if(aL.equals("Be")) number=number+4;
				if(aL.equals("B"))  number=number+5;
				if(aL.equals("C"))  number=number+6;
				if(aL.equals("N"))  number=number+7;
				if(aL.equals("O"))  number=number+8;
				if(aL.equals("F"))  number=number+9;
				if(aL.equals("Ne")) number=number+10;
				if(aL.equals("Na")) number=number+11;
				if(aL.equals("Mg")) number=number+12;
				if(aL.equals("Al")) number=number+13;
				if(aL.equals("Si")) number=number+14;
				if(aL.equals("P"))  number=number+15;
				if(aL.equals("S"))  number=number+16;
				if(aL.equals("Cl")) number=number+17;
				if(aL.equals("Ar")) number=number+18;
				if(aL.equals("K"))  number=number+19;
				if(aL.equals("Ca")) number=number+20;
				if(aL.equals("Sc")) number=number+21;
				if(aL.equals("Ti")) number=number+22;
				if(aL.equals("V"))  number=number+23;
				if(aL.equals("Cr")) number=number+24;
				if(aL.equals("Mn")) number=number+25;
				if(aL.equals("Fe")) number=number+26;
				if(aL.equals("Co")) number=number+27;
				if(aL.equals("Ni")) number=number+28;
				if(aL.equals("Cu")) number=number+29;
				if(aL.equals("Zn")) number=number+30;
				if(aL.equals("Ga")) number=number+31;
				if(aL.equals("Ge")) number=number+32;
				if(aL.equals("As")) number=number+33;
				if(aL.equals("Se")) number=number+34;
				if(aL.equals("Br")) number=number+35;
				if(aL.equals("Kr")) number=number+36;
				if(aL.equals("Rb")) number=number+37;
				if(aL.equals("Sr")) number=number+38;
				if(aL.equals("Y"))  number=number+39;
				if(aL.equals("Zr")) number=number+40;
				if(aL.equals("Nb")) number=number+41;
				if(aL.equals("Mo")) number=number+42;
				if(aL.equals("Tc")) number=number+43;
				if(aL.equals("Ru")) number=number+44;
				if(aL.equals("Rh")) number=number+45;
				if(aL.equals("Pd")) number=number+46;
				if(aL.equals("Ag")) number=number+47;
				if(aL.equals("Cd")) number=number+48;
				if(aL.equals("In")) number=number+49;
				if(aL.equals("Sn")) number=number+50;
				if(aL.equals("Sb")) number=number+51;
				if(aL.equals("Te")) number=number+52;
				if(aL.equals("I"))  number=number+53;
				if(aL.equals("Xe")) number=number+54;
				if(aL.equals("Cs")) number=number+55;
				if(aL.equals("Ba")) number=number+56;
				if(aL.equals("La")) number=number+57;
				if(aL.equals("Ce")) number=number+68;
				if(aL.equals("Pr")) number=number+69;
				if(aL.equals("Nd")) number=number+60;
				if(aL.equals("Pm")) number=number+61;
				if(aL.equals("Sm")) number=number+62;
				if(aL.equals("Eu")) number=number+63;
				if(aL.equals("Gd")) number=number+64;
				if(aL.equals("Tb")) number=number+65;
				if(aL.equals("Dy")) number=number+66;
				if(aL.equals("Ho")) number=number+67;
				if(aL.equals("Er")) number=number+68;
				if(aL.equals("Tm")) number=number+69;
				if(aL.equals("Yb")) number=number+70;
				if(aL.equals("Lu")) number=number+71;
				if(aL.equals("Hf")) number=number+72;
				if(aL.equals("Ta")) number=number+73;
				if(aL.equals("W"))  number=number+74;
				if(aL.equals("Re")) number=number+75;
				if(aL.equals("Os")) number=number+76;
				if(aL.equals("Ir")) number=number+77;
				if(aL.equals("Pt")) number=number+78;
				if(aL.equals("Au")) number=number+79;
				if(aL.equals("Hg")) number=number+80;
				if(aL.equals("Tl")) number=number+81;
				if(aL.equals("Pb")) number=number+82;
				if(aL.equals("Bi")) number=number+83;
				if(aL.equals("Po")) number=number+84;
				if(aL.equals("At")) number=number+85;
				if(aL.equals("Rn")) number=number+86;
			}
			br.close();
			return number;
		} catch(IOException ef) {
			System.out.println("IOException in XYZ file "+fileName);
			return 0;
		}
	}
	
	private void dimers() {
		
		Integer n=0;
		for(int i=0; i<numFragments-1; i++) {
			for(int j=i+1; j<numFragments; j++) {
				n++;
			}
		}
		newMEBFs=n;
		dimensionData[2][1]=newMEBFs;
		update();
		
		n=0;
		for(int i=0; i<numFragments-1; i++) {
			for(int j=i+1; j<numFragments; j++) {
				mebfFragments[n][0][0]=i;
				mebfFragments[n][1][0]=j;
				n++;
			}
		}

		if(numMEBFs>0) {
			for(int i=0; i<numMEBFs; i++) {
				Integer mebf = mebfIndex[i][0];
				mebf=i;
				Integer nmer = mebfSpecification[i][0];
				Integer spin = mebfSpecification[i][1];
				Integer chrg = mebfSpecification[i][2];
				Integer stat = mebfSpecification[i][3];
				mebfSpecification[i][0]=2;
				selectMEBFStates(mebf,nmer,spin,chrg,stat);	
			}
		}
		
		n=0;
		for(int i=0; i<numFragments-1; i++) {
			for(int j=i+1; j<numFragments; j++) {
				mebfFragments[n][0][0]=i;
				mebfFragments[n][1][0]=j;
				n++;
			}
		}

	}

	private void dimers_noA() {
		
		Integer n=0;
		for(int i=1; i<numFragments-1; i++) {
			for(int j=i+1; j<numFragments; j++) {
				n++;
			}
		}
		newMEBFs=n;
		dimensionData[2][1]=newMEBFs;
		update();
		
		n=0;
		for(int i=1; i<numFragments-1; i++) {
			for(int j=i+1; j<numFragments; j++) {
				mebfFragments[n][0][0]=i;
				mebfFragments[n][1][0]=j;
				n++;
			}
		}

		if(numMEBFs>0) {
			for(int i=0; i<numMEBFs; i++) {
				Integer mebf = mebfIndex[i][0];
				mebf=i;
				Integer nmer = mebfSpecification[i][0];
				Integer spin = mebfSpecification[i][1];
				Integer chrg = mebfSpecification[i][2];
				Integer stat = mebfSpecification[i][3];
				mebfSpecification[i][0]=2;
				selectMEBFStates(mebf,nmer,spin,chrg,stat);	
			}
		}
		
		n=0;
		for(int i=1; i<numFragments-1; i++) {
			for(int j=i+1; j<numFragments; j++) {
				mebfFragments[n][0][0]=i;
				mebfFragments[n][1][0]=j;
				n++;
			}
		}

	}
	private void update() {
		
		updateStatesList();
		updateFragmentList();
		updateFragmentDefinitions();
		readOutputFiles();
		updateFragmentEnergies();
		updateMEBFDefinitions();
		for(int i=0; i<numMEBFs; i++) {
			Integer mebf = mebfIndex[i][0];
			Integer nmer = mebfSpecification[mebf][0];
			Integer spin = mebfSpecification[mebf][1];
			Integer chrg = mebfSpecification[mebf][2];
			Integer stat = mebfSpecification[mebf][3];
			if(stat<=0) selectMEBFStates(mebf,nmer,spin,chrg,stat);
		}
		
		updateNOCIEnergies();
		updateHints();
		
		statesPanel.setPreferredSize(new Dimension(Short.MAX_VALUE,((numSets)*15+50)));
		statesPanel.setMinimumSize(new Dimension(Short.MAX_VALUE,((numSets)*15+50)));
		statesPanel.setMaximumSize(new Dimension(Short.MAX_VALUE,((numSets)*15+50)));
		
		fragmentsPanel.setPreferredSize(new Dimension(Short.MAX_VALUE,((numFragments)*15+50)));
		fragmentsPanel.setMinimumSize(new Dimension(Short.MAX_VALUE,((numFragments)*15+50)));
		fragmentsPanel.setMaximumSize(new Dimension(Short.MAX_VALUE,((numFragments)*15+50)));
		
		energiesPanel.setPreferredSize(new Dimension(Short.MAX_VALUE,((numEnergies)*15+50)));
		energiesPanel.setMinimumSize(new Dimension(Short.MAX_VALUE,((numEnergies)*15+50)));
		energiesPanel.setMaximumSize(new Dimension(Short.MAX_VALUE,((numEnergies)*15+50)));
		
		mebfsPanel.setPreferredSize(new Dimension(Short.MAX_VALUE,((numMEBFRows)*15+50)));
		mebfsPanel.setMinimumSize(new Dimension(Short.MAX_VALUE,((numMEBFRows)*15+50)));
		mebfsPanel.setMaximumSize(new Dimension(Short.MAX_VALUE,((numMEBFRows)*15+50)));
		
		nociResultsPanel.setPreferredSize(new Dimension(Short.MAX_VALUE,255));
		nociResultsPanel.setMinimumSize(new Dimension(Short.MAX_VALUE,255));
		nociResultsPanel.setMaximumSize(new Dimension(Short.MAX_VALUE,noci1*15+65));
		nociEnergiesPanel.setPreferredSize(new Dimension(noci0*80+40,noci1*15+50));
		nociEnergiesPanel.setMinimumSize(new Dimension(noci0*80+40,noci1*15+50));
		nociEnergiesPanel.setMaximumSize(new Dimension(noci0*80+40,noci1*15+50));
//		nociPlotPanel.setPreferredSize(new Dimension(300,250));
//		nociPlotPanel.setMinimumSize(new Dimension(300,250));
//		nociPlotPanel.setMaximumSize(new Dimension(300,250));
		nociPlotPanel.setPreferredSize(new Dimension(300,noci1*15+50));
		nociPlotPanel.setMinimumSize(new Dimension(300,noci1*15+50));
		nociPlotPanel.setMaximumSize(new Dimension(300,noci1*15+50));
		
		parametersPanel.revalidate();
		dimensionPanel.revalidate();
		numberPanel.revalidate();	
		threshPanel.revalidate();	
		statesPanel.revalidate();
		fragmentsPanel.revalidate();
		energiesPanel.revalidate();
		mebfsPanel.revalidate();
		nociEnergiesPanel.revalidate();
		nociPlotPanel.revalidate();
		nociResultsPanel.revalidate();
		
		parametersPanel.repaint();
		dimensionPanel.repaint();
		jobPanel.repaint();
		numberPanel.repaint();
		threshPanel.repaint();
		statesPanel.repaint();
		fragmentsPanel.repaint();
		energiesPanel.repaint();
		mebfsPanel.repaint();
		nociEnergiesPanel.repaint();
		nociPlotPanel.repaint();
		nociResultsPanel.repaint();
		
		container.revalidate();
		container.repaint();
		container.setVisible(true);

	}
	
	private void updateStatesList() {

		if(newSets>numSets) {
			System.out.println("New State List");
			for(int i=numSets; i<newSets; i++) {
		    	lenStateList[i]=0;
		    	spinStateList[i]=0;
		    	for(int j=0; j<15; j++) ndxStateList[i][j]=0;
			}
		}

		Integer[] count = new Integer[numSets];
		for(int i=0; i<numSets; i++) {
			count[i]=0;
			for(int j=0; j<numFragments; j++) {
				if(spinStateList[i]==1 && (dimFragments[i][3] % 2)==0) count[i]++;
				if(spinStateList[i]==2 && (dimFragments[i][3] % 2)==1) count[i]++;
			}
		}

		for(int i=0; i<numFragments; i++) {
			Boolean addSet = true;
			for(int j=0; j<numSets; j++) {
				if(spinStateList[j]==1 && (dimFragments[i][3] % 2)==0) addSet=false;
				if(spinStateList[j]==2 && (dimFragments[i][3] % 2)==1) addSet=false;
			}
			if(addSet) {
				for(int j=0; j<numSets; j++) {
					if(count[j]==0) {
				    	spinStateList[j]=(dimFragments[i][3] % 2)+1;
						lenStateList[j]=0;
						dimFragments[i][2]=j;
						addSet=false;
					}
				}
				if(addSet) {
					if(newSets==numSets) newSets++;
					spinStateList[numSets]=(dimFragments[i][3] % 2)+1;
					lenStateList[numSets]=0;
					dimFragments[i][2]=numSets;
				}
			}
		}

		if(numSets>numberStateDefinitions) {
			for(int i=numberStateDefinitions; i<numSets; i++) {
				statesTableModel.addRow(stateDefinitions[i]);
			}
			numberStateDefinitions=numSets;
		}

		if(newSets>numberStateDefinitions) {
			for(int i=numberStateDefinitions; i<newSets; i++) {
				statesTableModel.addRow(stateDefinitions[i]);
			}
		}

		if(newSets>0 && newSets<numberStateDefinitions) {
			for(int i=newSets; i<numberStateDefinitions; i++) {
				statesTableModel.removeRow(numberStateDefinitions-1-i);
			}
		}
		
		numSets=newSets;
		
		for(int i=0; i<numSets; i++) {
			Boolean redo = false;
			if(spinStateList[i]==1 || lenStateList[i]==0) {
				redo=false;
				for(int j=0; j<lenStateList[i]; j++) {
					if(stateNames[ndxStateList[i][j]].trim().equals("D0")) redo=true;
					if(stateNames[ndxStateList[i][j]].trim().equals("D1")) redo=true;
					if(stateNames[ndxStateList[i][j]].trim().equals("S+")) redo=true;
					if(stateNames[ndxStateList[i][j]].trim().equals("S-")) redo=true;
					if(stateNames[ndxStateList[i][j]].trim().equals("T+")) redo=true;
					if(stateNames[ndxStateList[i][j]].trim().equals("T-")) redo=true;
				}
				if(redo || lenStateList[i]==0) {
					lenStateList[i]=5;
					ndxStateList[i][0]=0;
					ndxStateList[i][1]=1;
					ndxStateList[i][2]=5;
					ndxStateList[i][3]=8;
					ndxStateList[i][4]=11;
				}
			}
			if(spinStateList[i]==2 || lenStateList[i]==0) {
				redo=false;
				for(int j=0; j<lenStateList[i]; j++) {
					if(stateNames[ndxStateList[i][j]].trim().equals("S0")) redo=true;
					if(stateNames[ndxStateList[i][j]].trim().equals("S1")) redo=true;
					if(stateNames[ndxStateList[i][j]].trim().equals("S2")) redo=true;
					if(stateNames[ndxStateList[i][j]].trim().equals("T1")) redo=true;
					if(stateNames[ndxStateList[i][j]].trim().equals("T2")) redo=true;
					if(stateNames[ndxStateList[i][j]].trim().equals("D+")) redo=true;
					if(stateNames[ndxStateList[i][j]].trim().equals("D-")) redo=true;
				}
				if(redo || lenStateList[i]==0) {
					lenStateList[i]=6;
					ndxStateList[i][0]=3;
					ndxStateList[i][1]=4;
					ndxStateList[i][2]=7;
					ndxStateList[i][3]=9;
					ndxStateList[i][4]=10;
					ndxStateList[i][5]=12;					
				}
			}
		}
	    String[] stateNames = new String[] {"S0","S1","S2","D0","D1","T1","T2","S+","D+","T+","S-","D-","T-","Q1"};

		for(int i=0; i<numSets; i++) {
			for(int j=0; j<18; j++) stateDefinitions[i][j]=" ";
			stateDefinitions[i][0]=(i+1);
			stateDefinitions[i][1]=lenStateList[i];
			stateDefinitions[i][2]=spinStateList[i];
			for(int j=0; j<lenStateList[i]; j++) stateDefinitions[i][j+3]=stateNames[ndxStateList[i][j]];
		}
		for(int i=0; i<numSets; i++) {
			for(int j=0; j<18; j++) {
				statesTableModel.setValueAt(stateDefinitions[i][j],i,j);
			}
		}
		numberStateDefinitions=numSets;	
	}
	
	private void updateFragmentList() {

		if(numFragments>numberFragmentDefinitions) {
			for(int i=numberFragmentDefinitions; i<numFragments; i++) {
				fragmentsTableModel.addRow(fragmentDefinitions[i]);
			}
			numberFragmentDefinitions=numFragments;
		}

		if(newFragments>numberFragmentDefinitions) {
			for(int i=numberFragmentDefinitions; i<newFragments; i++) {
				fragmentsTableModel.addRow(fragmentDefinitions[i]);
			}
		}

		if(newFragments>0 && newFragments<numberFragmentDefinitions) {
			for(int i=newFragments; i<numberFragmentDefinitions; i++) {
				fragmentsTableModel.removeRow(numberFragmentDefinitions-1-i);
			}
		}

		if(newFragments>numFragments) {			
			for(int i=numFragments; i<newFragments; i++) {
		    	namFragments[i]="";
		    	for(int j=0; j<6; j++) movFragments[i][j]=0.0;
		    	for(int j=0; j<11; j++) dimFragments[i][j]=0;
		    	energiesDFT[i]=0.0;
		    	energiesSCF[i]=0.0;
		    	for(int j=0; j<7; j++) {
		    		energiesCASSCF[i][j]=0.0;
		    		energiesCASPT2[i][j]=0.0;
		    	}
		    	dimFragments[i][2]=0;
		    	dimFragments[i][4]=numCASe;
		    	dimFragments[i][5]=numCASo;
			}
			if(numFragments>0) {		
				for(int i=numFragments; i<newFragments; i++) {
					namFragments[i]=namFragments[numFragments-1];
					dimFragments[i][1]=dimFragments[numFragments-1][1];
					dimFragments[i][2]=dimFragments[numFragments-1][2];
					dimFragments[i][3]=dimFragments[numFragments-1][3];
					dimFragments[i][4]=dimFragments[numFragments-1][4];
					dimFragments[i][5]=dimFragments[numFragments-1][5];
				}
			}
		}
		
		numFragments=newFragments;

		for(int i=0; i<numFragments; i++) {
			fragmentDefinitions[i][0]=fragmentNames[i];
			if(fragmentDefinitions[i][1]==null) fragmentDefinitions[i][1]=fragmentDefinitions[dimFragments[i][0]][0];
			fragmentDefinitions[i][2]=namFragments[i];
			fragmentDefinitions[i][3]=dimFragments[i][1];
			fragmentDefinitions[i][4]=dimFragments[i][2]+1;
			fragmentDefinitions[i][5]=dimFragments[i][3];
			fragmentDefinitions[i][6]=dimFragments[i][4];
			fragmentDefinitions[i][7]=dimFragments[i][5];
			fragmentDefinitions[i][8]=movFragments[i][0];
			fragmentDefinitions[i][9]=movFragments[i][1];
			fragmentDefinitions[i][10]=movFragments[i][2];
			fragmentDefinitions[i][11]=movFragments[i][3];
			fragmentDefinitions[i][12]=movFragments[i][4];
			fragmentDefinitions[i][13]=movFragments[i][5];
			fragmentDefinitions[i][14]=dimFragments[i][10];
		}

		for(int i=0; i<numFragments; i++) {
			for(int j=0; j<15; j++) {
				fragmentsTableModel.setValueAt(fragmentDefinitions[i][j],i,j);
			}
		}
		numberFragmentDefinitions=numFragments;
		
		// Update the Fragment Energies table headers

		numFragmentStates=0;
		for(int i=0; i<numStateLabels; i++) {
			stateLabels[i]=" ";
			stateLabelIndex[i]=i;
		}
		for(int i=0; i<numFragments; i++) {
			Integer stateIndex=dimFragments[i][2];
				if(lenStateList[stateIndex]>0) {
				for(int j=0; j<lenStateList[stateIndex]; j++) {
					Boolean found = false;
					for(int k=0; k<numFragmentStates; k++) {
						if(stateNames[ndxStateList[stateIndex][j]].trim().equals(stateLabels[(k+1)].trim())) found=true;
					}
					if(!found) {
						numFragmentStates++;
						stateLabels[numFragmentStates]=stateNames[ndxStateList[stateIndex][j]].trim();
						stateLabelIndex[ndxStateList[stateIndex][j]]=numFragmentStates;
					}
				}
			}
		}
		for(int i=0; i<12; i++) {
			energiesTable.getColumnModel().getColumn(i).setHeaderValue(stateLabels[i]);
		}
	}
	
	private void updateFragmentDefinitions() {

		String nameA, nameB, nameP, nameF;
		
		for(int i=0; i<numEnergies; i++) {
			for(int j=0; j<8; j++) {
				stateEnergies[i][j]=null;
				energiesTableModel.setValueAt(stateEnergies[i][j],i,j);
			}
		}
		
		for(int i=0; i<numFragments; i++) dimFragments[i][0]=i;
		for(int i=0; i<numFragments-1; i++) {
			for(int j=i+1; j<numFragments; j++) {
				if(namFragments[i].trim().equals(namFragments[j].trim())) {
					dimFragments[j][0]=dimFragments[i][0];	
				}
			}
		}
		
		for(int i=0; i<numFragments; i++) {
			if(fragmentDefinitions[i][1]==null) fragmentDefinitions[i][1]=fragmentDefinitions[dimFragments[i][0]][0];
			fragmentDefinitions[i][2]=namFragments[i];
			fragmentDefinitions[i][3]=dimFragments[i][1];
			fragmentDefinitions[i][4]=dimFragments[i][2]+1;
			fragmentDefinitions[i][5]=dimFragments[i][3];
			fragmentDefinitions[i][6]=dimFragments[i][4];
			fragmentDefinitions[i][7]=dimFragments[i][5];
			fragmentDefinitions[i][8]=movFragments[i][0];
			fragmentDefinitions[i][9]=movFragments[i][1];
			fragmentDefinitions[i][10]=movFragments[i][2];
			fragmentDefinitions[i][11]=movFragments[i][3];
			fragmentDefinitions[i][12]=movFragments[i][4];
			fragmentDefinitions[i][13]=movFragments[i][5];
			fragmentDefinitions[i][14]=dimFragments[i][10];
		}
		
		numEnergies=0;
		Integer index=0;
		for(int i=0; i<numFragments; i++ ) {
			if(dimFragments[i][0]==i) {
				numEnergies=numEnergies+4;
				for(int j=0; j<4; j++) {
					stateEnergies[index][0]=(String) fragmentNames[i]+":"+energyNames[j];
					energyFragment[index][0]=i;
					energyFragment[index][1]=j;
					index++;
				}
			} else {
				numEnergies=numEnergies+1;
				for(int j=2; j<3; j++) {
					stateEnergies[index][0]=(String) fragmentNames[i]+":"+energyNames[j];
					energyFragment[index][0]=i;
					energyFragment[index][1]=j;
					index++;
				}
			}
		}
		
		if(numEnergies>numberStateEnergies) {
			for(int i=numberStateEnergies; i<numEnergies; i++) {
				energiesTableModel.addRow(stateEnergies[i]);
			}
		}
		
		if(numEnergies<numberStateEnergies) {
			for(int i=numEnergies; i<numberStateEnergies; i++) {
				energiesTableModel.removeRow(numberStateEnergies-1-i);
			}
		}

		for(int i=0; i<numEnergies; i++) {
			energiesTableModel.setValueAt(stateEnergies[i][0],i,0);
		}

		for(int i=0; i<numFragments; i++) {
			for(int j=0; j<6; j++) RandT[j]=movFragments[i][j];
			nameP=projectName.trim();
			nameF=namFragments[i].trim();
			nameB=(String) fragmentDefinitions[dimFragments[i][0]][0];
			nameA=fragmentNames[i].trim();
			if(nameB.trim().equals(nameA.trim())) {
				nameB=" ";
			} else {
//				Why setting fragment name to project name???
//				nameF=nameP.trim();
			}
			fragment.initialize(nameP,nameF,nameA,nameB,dimFragments[i][3],RandT);

		}
		numberStateEnergies=numEnergies;
	}

	private void updateFragmentEnergies() {
		Integer ifrag;
		Integer itype;
		Integer istat;
		Double energy;
		Integer maxCol=1;

		for(int i=0; i<numEnergies; i++) {
			for(int j=0; j<12; j++) stateEnergies[i][j]=" ";
		}
		
		for(int i=0; i<numFragments-1; i++) {
				for(int j=i+1; j<numFragments; j++) {
					if(!fragmentDefinitions[j][0].equals(fragmentDefinitions[j][1])) {
						if(fragmentDefinitions[j][1].equals(fragmentDefinitions[i][0])) {
							istat=dimFragments[i][2];
							for(int k=0; k<lenStateList[istat]; k++) {
							energiesCASSCF[j][k]=energiesCASSCF[i][k];
							energiesCASPT2[j][k]=energiesCASPT2[i][k];
						}
					}
				}
			}
		}
		
		for(int i=0; i<numEnergies; i++) {
			// fragment = energyFragment[i][0]
			// energy type = energyFragment[i][1] where 0=DFT 1=SCF 2=CASSCF 3=CASPT2
			ifrag=energyFragment[i][0];
			itype=energyFragment[i][1];
			istat=dimFragments[ifrag][2];
			if(itype==0) {
				energy=energiesDFT[ifrag];
				if(energy!=0.0) stateEnergies[i][0]=energy;
			}
			if(itype==1) {
				energy=energiesSCF[ifrag];
				if(energy!=0.0) stateEnergies[i][0]=energy;
			}
			if(itype==2 && lenStateList[istat]>0) {
				for(int j=0; j<lenStateList[istat]; j++) {
					energy=energiesCASSCF[ifrag][j];					
					if(energy!=0.0) stateEnergies[i][stateLabelIndex[ndxStateList[istat][j]]-1]=energy;
					maxCol=Math.max(maxCol,j);
				}
			}
			if(itype==3) {
				for(int j=0; j<lenStateList[istat]; j++) {
					energy=energiesCASPT2[ifrag][j];
					if(energy!=0.0) stateEnergies[i][stateLabelIndex[ndxStateList[istat][j]]-1]=energy;					
					maxCol=Math.max(maxCol,j);
				}
			}
		}
		for(int i=0; i<numEnergies; i++) {
			for(int j=0; j<=maxCol; j++) energiesTableModel.setValueAt(stateEnergies[i][j],i,j+1);
		}
		energiesTable.repaint();
	}

	public void updateNOCIEnergies() {

		Integer ndx=0;
		
		if(noci0<2*numMEBFs+1) {
			for(int i=noci0; i<2*numMEBFs+1; i=i+2) {
				ndx=(i-1)/2;
				nociEnergiesTableModel.addColumn("E_MEBF("+mebfName[ndx]+")");
				nociEnergiesTableModel.addColumn("E_NOCI("+mebfName[ndx]+")");
			}
			noci0=2*numMEBFs+1;
		}
		if(noci0>2*numMEBFs+1) {
			noci0=2*numMEBFs+1;
			nociEnergiesTableModel.setColumnCount(noci0);
		}

		for(int i=1; i<2*numMEBFs+1; i=i+2) {
			ndx=(i-1)/2;
			
			nociEnergiesTable.getColumnModel().getColumn(i).setHeaderValue("E_MEBF("+mebfName[ndx]+")");
			nociEnergiesTable.getColumnModel().getColumn(i+1).setHeaderValue("E_NOCI("+mebfName[ndx]+")");
		}
		
		Integer maxMEBFStates=0;
		
		for(int i=0; i<numMEBFs; i++) {
			if(maxMEBFStates<mebfSpecification[i][3]) maxMEBFStates=mebfSpecification[i][3];
		}
		
		if(maxMEBFStates>0) {	
			if(noci1<maxMEBFStates) {
				for(int i=noci1; i<maxMEBFStates; i++) {
					nociEnergiesTableModel.addRow(nociEnergies[i-1]);
				}
				noci1=maxMEBFStates;
			}
			
			if(noci1>maxMEBFStates) {
				noci1=maxMEBFStates;
				nociEnergiesTableModel.setRowCount(noci1);
			}
			
			if(noci0>0) {
				for(int i=0; i<noci0; i++) {
					for(int j=0; j<noci1; j++) {
						nociEnergies[i][j]=" ";
						nociEnergiesTableModel.setValueAt(nociEnergies[i][j],j,i);
					}
				}
			}
	
			for(int j=0; j<numMEBFs; j++) {
				for(int i=0; i<mebfSpecification[j][3]; i++) {
					nociEnergies[i][j]=" ";
				}
			}
			
			nociEnergiesTable.getColumnModel().getColumn(0).setMaxWidth(40);
			for(int j=1; j<2*numMEBFs+1; j++) nociEnergiesTable.getColumnModel().getColumn(j).setMaxWidth(100);

			energyDiagram.clear(false);

			for(int i=0; i<noci1; i++) nociEnergiesTableModel.setValueAt(" ",i,0);
			for(int i=0; i<maxMEBFStates; i++) nociEnergiesTableModel.setValueAt("E"+(i+1),i,0);
			
			Integer nums=0;
			
			for(int i=0; i<numMEBFs; i++) {
				if(read_GronOR_arx(i)) {
					if(numEnergiesGronOR[i][0]>0) {
						for(int j=0; j<numEnergiesGronOR[i][0]; j++) {
							nociEnergies[j][2*i]=energiesGronOR[j][i][0];
							nociEnergiesTableModel.setValueAt(nociEnergies[j][2*i],j,2*i+1);
						}
					}
					if(numEnergiesGronOR[i][1]>0) {
						for(int j=0; j<numEnergiesGronOR[i][1]; j++) {
							nociEnergies[j][2*i+1]=energiesGronOR[j][i][1];
							nociEnergiesTableModel.setValueAt(nociEnergies[j][2*i+1],j,2*i+2);
						}
					}
					nums=i+1;
				}
			}
			
			for(int i=0; i<maxMEBFStates; i++) nociEnergiesTableModel.setValueAt("E"+(i+1),i,0);
				
			if(nums>0) {				
				Double Emin=0.0;
				Double Emax=0.0;
				for(int i=0; i<nums; i++) {
					for(int j=0; j<numEnergiesGronOR[i][0]; j++) {
						Emin=Math.min(Emin, energiesGronOR[j][i][0]);
					}
					for(int j=0; j<numEnergiesGronOR[i][1]; j++) {
						Emin=Math.min(Emin, energiesGronOR[j][i][1]);
					}
				}
				Emax=Emin;
				for(int i=0; i<nums; i++) {
					for(int j=0; j<numEnergiesGronOR[i][0]; j++) {
						Emax=Math.max(Emax, energiesGronOR[j][i][0]);
					}
					for(int j=0; j<numEnergiesGronOR[i][1]; j++) {
						Emax=Math.max(Emax, energiesGronOR[j][i][1]);
					}
				}
				for(int i=0; i<nums; i++) {
					for(int j=0; j<numEnergiesGronOR[i][0]; j++) {
						energiesRelGronOR[j][i][0]=energiesGronOR[j][i][0]-Emin;
					}
					for(int j=0; j<numEnergiesGronOR[i][1]; j++) {
						energiesRelGronOR[j][i][1]=energiesGronOR[j][i][1]-Emin;
					}
				}
				
				energyDiagram.setForeground(Color.red);
				
				Double Edelta = Emax-Emin;
				energyDiagram.setYRange(-0.1*Edelta,1.1*Edelta);
				energyDiagram.setXRange(0.0,Double.valueOf(4*numMEBFs-1));
				for(int i=0; i<nums; i++) {
					for(int j=0; j<numEnergiesGronOR[i][0]; j++) {
						energyDiagram.addPoint(0,Double.valueOf(4*i),energiesRelGronOR[j][i][0],false);
						energyDiagram.addPoint(0,Double.valueOf(4*i+1),energiesRelGronOR[j][i][0],true);
					}
					for(int j=0; j<numEnergiesGronOR[i][1]; j++) {
						energyDiagram.addPoint(0,Double.valueOf(4*i+2),energiesRelGronOR[j][i][1],false);
						energyDiagram.addPoint(0,Double.valueOf(4*i+3),energiesRelGronOR[j][i][1],true);
					}
				}
				energyDiagram.setForeground(Color.black);
				for(int i=0; i<nums; i++) {
					for(int j=0; j<numEnergiesGronOR[i][1]; j++) {
						for(int k=0; k<numEnergiesGronOR[i][0]; k++) {
							if(Math.abs(coefGronOR[k][j][i])>0.75) {
								energyDiagram.addPoint(0,Double.valueOf(4*i+1),energiesRelGronOR[k][i][0],false);
								energyDiagram.addPoint(0,Double.valueOf(4*i+2),energiesRelGronOR[j][i][1],true);
							}
						}
					}
				}
			}
		}
//		energyDiagram.revalidate();
//		energyDiagram.repaint();
	}
	
	private void updateHints() {
		Boolean fragT[] = new Boolean[numFragments];
		Integer nfrag = 0;
		if(numSets>0) setTableCellColor(dimensionTable,0,1,Color.white);
		if(numFragments>0) {
			setTableCellColor(dimensionTable,1,1,Color.white);
			setTableCellColor(fragmentsTable,0,2,Color.white);
			clearTableColor(fragmentsTable,0,2,Color.white);
			for(int i=0; i<numFragments; i++) {
				fragT[i]=false;
				if(movFragments[i][0]!=0.0 || movFragments[i][1]!=0.0 || movFragments[i][2]!=0.0) fragT[i]=true;
				setTableCellColor(fragmentsTable,i,8,Color.white);
				setTableCellColor(fragmentsTable,i,9,Color.white);
				setTableCellColor(fragmentsTable,i,10,Color.white);
				if(!fragT[i]) nfrag++;
			}
		}
		if(numMEBFs>0) setTableCellColor(dimensionTable,2,1,Color.white);
		clearTableColor(dimensionTable,2,1,Color.white);
		if(numSets<=0) {
			setTableCellColor(dimensionTable,0,1,Color.red);
		} else if(numFragments<=0) {
			setTableCellColor(dimensionTable,1,1,Color.red);
		} else if(fragmentDefinitions[0][2]=="") {
			setTableCellColor(fragmentsTable,0,2,Color.red);
		} else if(nfrag>1) {
			for(int i=1; i<numFragments; i++) {
				if(!fragT[i]) {
					setTableCellColor(fragmentsTable,i,8,Color.red);
					setTableCellColor(fragmentsTable,i,9,Color.red);
					setTableCellColor(fragmentsTable,i,10,Color.red);
				}
			}
		} else if(numMEBFs<=0){
			setTableCellColor(dimensionTable,2,1,Color.red);
		}
	}
	
	private void setTableCellColor(JTable table, Integer row, Integer column, Color color) {
		table.setCellSelectionEnabled(true);
		TableCellRenderer renderer = table.getCellRenderer(row, column);
		Object value = table.getModel().getValueAt(row, column);
		Boolean selected = table.getSelectionModel().isSelectedIndex(row);
		Boolean focus = false;
		selected = false;
		Component component= (JLabel) renderer.getTableCellRendererComponent(table, value, selected, focus, row, column);
		table.setRowSelectionInterval(row, row);
		table.setColumnSelectionInterval(column, column);
		table.setSelectionBackground(color);
	}
	
	private void clearTableColor(JTable table, Integer row, Integer column, Color color) {
		table.setCellSelectionEnabled(true);
		TableCellRenderer renderer = table.getCellRenderer(row, column);
		Object value = table.getModel().getValueAt(row, column);
		Boolean selected = table.getSelectionModel().isSelectedIndex(row);
		Boolean focus = false;
		selected = false;
		Component component= (JLabel) renderer.getTableCellRendererComponent(table, value, selected, focus, row, column);
		component.setBackground(color);		
	}
	
	private void writeInputFiles() {

		String nameA, nameP, nameF, nameS;
		Integer numCASe=0, numCASo=0;
		Boolean withCASPT2=false;
		Boolean withAlter=true;
		
		Integer stateIndex = 0;

		Integer numRunDFT = -1;
		Integer numRunMolcas = -1;

		for(int i=0; i<numEnergies; i++) {
			if(dimFragments[energyFragment[i][0]][0]>numRunDFT) numRunDFT = dimFragments[energyFragment[i][0]][0];
			if(energyFragment[i][0]>numRunMolcas) numRunMolcas = energyFragment[i][0];
		}

		for(int i=0; i<=numRunDFT; i++) {
			for(int j=0; j<6; j++) RandT[j]=movFragments[i][j];
			nameP=projectName.trim();
			nameA= (String) fragmentDefinitions[dimFragments[i][0]][0];
			Integer mult=1;
			if(stateNames[ndxStateList[dimFragments[i][2]][0]].startsWith("D")) mult=2;
			if(stateNames[ndxStateList[dimFragments[i][2]][0]].startsWith("T")) mult=3;
//			nameF=projectName.trim()+nameA.trim();
			int last=namFragments[i].indexOf("_");
			nameF=namFragments[i].substring(0,last)+nameA.trim();
			nameP=projectName.trim()+nameA.trim();
			fragment.write_NWChem_DFT(nameF,nameP,mult,numRanks,RandT,memory,account,jobName,timeLimit);
		}

		for(int i=0; i<=numRunMolcas; i++) {
			for(int j=0; j<6; j++) RandT[j]=movFragments[i][j];
			nameP=projectName.trim();
			nameA= (String) fragmentDefinitions[i][0];
			stateIndex=dimFragments[i][2];
			withCASPT2=fragmentDefinitions[i][0].equals(fragmentDefinitions[i][1]);

			Integer mult=1;
			if(stateNames[ndxStateList[stateIndex][0]].startsWith("D")) mult=2;
			if(stateNames[ndxStateList[stateIndex][0]].startsWith("T")) mult=3;
			if(stateNames[ndxStateList[stateIndex][0]].startsWith("q")) mult=4;
			if(stateNames[ndxStateList[stateIndex][0]].startsWith("Q")) mult=5;
			int last=namFragments[i].indexOf("_");
			nameF=namFragments[i].substring(0,last)+nameA.trim();
			nameP=projectName.trim()+nameA.trim();
			fragment.write_Molcas_Int(nameF,nameP,basisSets[basisSet],contracts[contract],cholesky);
			fragment.write_Molcas_SCF(nameF,nameP,mult);

			if(lenStateList[stateIndex]>0) {
				withAlter=true;
				for(int j=0; j<lenStateList[stateIndex]; j++) {
					nameS=stateNames[ndxStateList[stateIndex][j]];
					numCASe=dimFragments[i][4];
					numCASo=dimFragments[i][5];
					fragment.write_Molcas_CASSCF(nameF,nameP,nameS,withCASPT2,numCASe,numCASo,withAlter,ipea);
					withAlter=false;
				}		
			}
		}

		for(int i=0; i<numFragments; i++) {
			nameA=(String) fragmentDefinitions[i][0];
			String nameB=(String) fragmentDefinitions[i][1];
			int last=namFragments[i].indexOf("_");
			nameF=namFragments[i].substring(0,last);
			nameP=projectName.trim();
			if(fragmentDefinitions[i][0].equals(fragmentDefinitions[i][1])) {
				fragment.write_Run_Script_Fragments(nameF,nameP,nameA,i,dimFragments[i][2],lenStateList,ndxStateList,numRanks,memory,account,jobName,timeLimit);
			} else {
				fragment.write_rotharm_input(nameP,nameA,nameB,dimFragments[i][2],lenStateList,ndxStateList);
				fragment.write_Rotate_Script_Fragments(nameF,nameP,nameA,i,dimFragments[i][2],lenStateList,ndxStateList,fragmentDefinitions[i][0],fragmentDefinitions[i][1],fragmentDefinitions[i][2],numRanks,memory,account,jobName,timeLimit);
			}
		}
		write_Molcas_MEBF_files();
		write_GronOR_NOCI();
		write_Run_All_Script();
	}
	
	private void write_Run_All_Script() {
		String fileName = projectName.trim()+"_all.run";
		try {
			PrintfWriter allFile = new PrintfWriter(new FileWriter(fileName));
			for(int i=0; i<numFragments; i++) {
				allFile.println("./"+projectName.trim()+fragmentNames[i].trim()+".run");
			}
			for(int i=0; i<numMEBFs; i++) {
				allFile.println("./"+projectName.trim()+mebfName[i].trim()+"_Molcas.run");
			}
			for(int i=0; i<numMEBFs; i++) {
				allFile.println("./"+projectName.trim()+mebfName[i].trim()+"_GronOR.run");
			}
			allFile.close();
		} catch(IOException ei) {
		}
	}
	
	private void readOutputFiles() {

		String nameA, nameB, nameP, nameF;
		
		Double energy = 0.0;
		Integer numAlt = 0;

		for(int i=0; i<numFragments; i++) {
			for(int j=0; j<6; j++) RandT[j]=movFragments[i][j];
			nameP=projectName.trim();
			nameF=namFragments[i].trim();
			nameA=" ";
			nameB=(String) fragmentDefinitions[dimFragments[i][0]][1];

			fragment.initialize(nameP, nameF, nameA, nameB, dimFragments[i][3], RandT);
			
			dimFragments[i][6]=0;
			dimFragments[i][7]=0;
			dimFragments[i][8]=0;
			dimFragments[i][9]=0;
			
			if(fragment.NWChem_Converged(i)) {
				energy=fragment.NWChem_DFT(i);
				energiesDFT[i]=energy;
				dimFragments[i][6]=1;
			} 
			
			nameA= (String) fragmentDefinitions[dimFragments[i][0]][0];
			numAlt=fragment.read_Alt(nameP,nameA);
			
			if(fragment.Molcas_SCF_Converged(i,dimFragments[i][4])) {
				energy=fragment.Molcas_SCF(i,dimFragments[i][4]);
				dimFragments[i][10]=numAlt;
				energiesSCF[i]=energy;
				dimFragments[i][7]=1;
			}
			
			for(int j=0; j<lenStateList[dimFragments[i][2]]; j++) {
				Integer stateIndex = ndxStateList[dimFragments[i][2]][j];
				Integer conv = fragment.Molcas_CASSCF_Converged(i,stateIndex);
				energiesCASSCF[i][j]=0.0;
				energiesCASPT2[i][j]=0.0;
				if(conv>=1) {
					energy=fragment.Molcas_CASSCF(i,stateIndex);
					energiesCASSCF[i][j]=energy;
					dimFragments[i][8]=j+1;
					if(conv==2) {
						energy=fragment.Molcas_CASPT2(i,stateIndex);
						energiesCASPT2[i][j]=energy;
						dimFragments[i][9]=j+1;
					}
				}
			}
		}		
//		read_GronOR_arx();
	}

	private void updateMEBFDefinitions() {

		Integer index = 0;
		// Add spin 1 1-mers as additional MEBFs
		if(newMEBFs>numMEBFs) {
			for(int i=numMEBFs; i<newMEBFs; i++) {
				mebfSpecification[i][0]=1;		// n-mer (i.e. number of fragments, default monomer)
				mebfSpecification[i][1]=1;		// spin  (default spin 1)
				mebfSpecification[i][2]=0;		// charge  (default charge 0)
				mebfSpecification[i][3]=0;		// number of states included (default one)
				mebfSpecification[i][4]=1;		// fragment index (default first)
				mebfName[i]=fragmentNames[0];	// mebf name (default name of first fragment)
				for(int j=0; j<maxMer; j++) {
					mebfFragments[i][j][0]=j;		// fragment index
					for(int k=1; k<maxMEBFstates; k++) {
						mebfFragments[i][j][k]=0;
					}
				}
				index++;
				if(index+1>numFragments) index=numFragments-1;
			}
		}

		//Count the number of required rows in MEBF table
		Integer oldRows = 0;
		for(int i=0; i<numMEBFs; i++) {
			oldRows=oldRows+mebfSpecification[i][0];
			mebfName[i]="";
			for(int j=0; j<mebfSpecification[i][0]; j++) {
				mebfName[i]=mebfName[i]+fragmentNames[mebfFragments[i][j][0]];
			}
		}
		
		//Count the number of required rows in MEBF table
		Integer numRows = 0;
		for(int i=0; i<newMEBFs; i++) {
			numRows=numRows+mebfSpecification[i][0];
			mebfName[i]="";
			for(int j=0; j<mebfSpecification[i][0]; j++) {
				mebfName[i]=mebfName[i]+fragmentNames[mebfFragments[i][j][0]];
			}
		}

		mebfsTableModel.setRowCount(numRows);

		numMEBFs=newMEBFs;

		for(int i=0; i<numRows; i++) {
			for(int k=0; k<maxMEBFstates+6; k++) {
				mebfDefinitions[i][k]=" ";
			}
		}
		
		
		if(numMEBFs>1) {
			for(int i=0; i<numMEBFs; i++) {
				index=1;
				for(int j=i+1; j<numMEBFs; j++ ) {
					if(mebfName[i].trim().equals(mebfName[j].trim())) {
						index++;
						mebfName[j]=mebfName[i].trim()+index;
					}
				}
				if(index>1) mebfName[i]=mebfName[i].trim()+"1";
			}
		}

		//Fill the mebfDefinitions table
		index=0;
		for(int i=0; i<numMEBFs; i++) {
			mebfDefinitions[index][0]=mebfName[i];					// mebf name
			mebfDefinitions[index][1]=mebfSpecification[i][0];		// n-mer
			mebfDefinitions[index][2]=mebfSpecification[i][1];		// spin
			mebfDefinitions[index][3]=mebfSpecification[i][2];		// charge
			mebfDefinitions[index][4]=mebfSpecification[i][3];		// number of states included
			// loop of number of fragments
			for(int j=0; j<mebfSpecification[i][0]; j++) {
				mebfDefinitions[index][5]=fragmentNames[mebfFragments[i][j][0]];
				// loop over number of states defined
				for(int k=0; k<mebfSpecification[i][3]; k++) {
					mebfDefinitions[index][6+k]=stateNames[mebfFragments[i][j][k+1]];
				}
				mebfIndex[index][0]=i;
				mebfIndex[index][1]=j;
				index++;
			}
		}
		//Fill the Table Model
		for(int i=0; i<numRows; i++) {
			for(int j=0; j<maxMEBFstates+6; j++) {
				mebfsTableModel.setValueAt(mebfDefinitions[i][j],i,j);
			}
		}

		numMEBFRows=numRows;
	}

	
	
	private void initializeWindow() {
		
		container = getContentPane();
		Box baseBox = Box.createVerticalBox();
		
		container.add(baseBox);
		 
// Parameter Panel
		
		parametersPanel = new JPanel();
		parametersPanel.setLayout(new BoxLayout(parametersPanel,BoxLayout.X_AXIS));
		parametersPanel.setPreferredSize(new Dimension(Short.MAX_VALUE,80));
		parametersPanel.setMinimumSize(new Dimension(Short.MAX_VALUE,80));
		parametersPanel.setMaximumSize(new Dimension(Short.MAX_VALUE,80));
		TitledBorder parametersBorder = new TitledBorder(new LineBorder(Color.black),"General Parameters");
		parametersBorder.setTitleColor(Color.black);
		parametersPanel.setBorder(parametersBorder);

		JPanel optionsPanel = new JPanel();
		optionsPanel.setLayout(new BoxLayout(optionsPanel,BoxLayout.Y_AXIS));
		optionsPanel.setPreferredSize(new Dimension(115,70));
		optionsPanel.setMinimumSize(new Dimension(115,70));
		optionsPanel.setMaximumSize(new Dimension(115,70));
		LineBorder optionsBorder = new LineBorder(Color.black);
		optionsPanel.setBorder(optionsBorder);
		JButton dimersButton = new JButton("Dimers");
		dimersButton.setPreferredSize(new Dimension(110,20));
		dimersButton.setMinimumSize(new Dimension(110,20));
		dimersButton.setMaximumSize(new Dimension(110,20));
		dimersButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				dimers();
				update();
				System.out.println("Constructed dimers");
			}
		});
		JButton dimersNoAButton = new JButton("Dimers -A");
		dimersNoAButton.setPreferredSize(new Dimension(110,20));
		dimersNoAButton.setMinimumSize(new Dimension(110,20));
		dimersNoAButton.setMaximumSize(new Dimension(110,20));
		dimersNoAButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				dimers_noA();
				update();
				System.out.println("Constructed dimers excluding A");
			}
		});
		
		dimensionPanel = new JPanel();
		dimensionPanel.setLayout(new BoxLayout(dimensionPanel,BoxLayout.X_AXIS));
		dimensionPanel.setPreferredSize(new Dimension(120,70));
		dimensionPanel.setMinimumSize(new Dimension(120,70));
		dimensionPanel.setMaximumSize(new Dimension(120,70));
		LineBorder dimensionBorder = new LineBorder(Color.black);
		dimensionPanel.setBorder(dimensionBorder);
		dimensionData = new Object[][] {
			{"State Sets",numSets},
			{"Fragments",numFragments},
			{"MEBFs",numMEBFs}
		};
		String[] dimensionColumns = new String[] {" "," "};
		dimensionTable = new JTable(dimensionData,dimensionColumns);
		dimensionTable.setCellSelectionEnabled(true);
		dimensionTable.getColumnModel().getColumn(0).setMaxWidth(80);
		dimensionTable.getColumnModel().getColumn(1).setMaxWidth(30);

		dimensionTable.addMouseListener(new MouseAdapter() {
			public void mouseClicked(java.awt.event.MouseEvent e) {
				Integer row=dimensionTable.getSelectedRow();
				Integer col=dimensionTable.getSelectedColumn();
				JFrame jf = new JFrame();
				String value;
				if(row==0 && col==1) {
					try {
						value = JOptionPane.showInputDialog(jf,"Enter new number of state definitions");
						if(value.length()>0) newSets=Integer.valueOf(value);
					} catch(NullPointerException e1) {
						newSets=numSets;
					}
//					update();
				}
				if(row==1 && col==1) {
					try {
						value = JOptionPane.showInputDialog(jf,"Enter new number of fragments");
						if(value.length()>0) newFragments=Integer.valueOf(value);
					} catch(NullPointerException e1) {
						newFragments=numFragments;
					}
//					update();
				}
				if(row==2 && col==1) {
					try {
						value = JOptionPane.showInputDialog(jf,"Enter new number of MEBFs");
						if(value.length()>0) newMEBFs=Integer.valueOf(value);
					} catch(NullPointerException e1) {
						newMEBFs=numMEBFs;
					}
//					update();
				}
				dimensionData[0][1]=newSets;
				dimensionData[1][1]=newFragments;
				dimensionData[2][1]=newMEBFs;
				update();

			}
		});
	
		
		numberPanel = new JPanel();
		numberPanel.setLayout(new BoxLayout(numberPanel,BoxLayout.X_AXIS));
		numberPanel.setPreferredSize(new Dimension(120,70));
		numberPanel.setMinimumSize(new Dimension(120,70));
		numberPanel.setMaximumSize(new Dimension(120,70));
		LineBorder numberBorder = new LineBorder(Color.black);
		numberPanel.setBorder(numberBorder);
		Object[][] numberData = new Object[][] {
			{"Ranks", numRanks},
			{"Threads", numThreads},
			{"Memory", memory}
		};
		String[] numberColumns = new String[] {" "," "};
		numberTable = new JTable(numberData,numberColumns);
		numberTable.setCellSelectionEnabled(true);
		numberTable.getColumnModel().getColumn(0).setMaxWidth(60);
		numberTable.getColumnModel().getColumn(1).setMaxWidth(50);
		ListSelectionModel numberSelectionModel = numberTable.getSelectionModel();
		numberSelectionModel.addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent e) {
				JFrame jf = new JFrame();
				String value;
				if(numberTable.getSelectedRow()==0) {
					try {
					value = JOptionPane.showInputDialog(jf,"Enter new number of ranks");
					if(value.length()>0) numRanks=Integer.valueOf(value);
					} catch(NullPointerException e1) {
					}
				}
				if(numberTable.getSelectedRow()==1) {
					try {
					value = JOptionPane.showInputDialog(jf,"Enter new number of threads");
					if(value.length()>0) numThreads=Integer.valueOf(value);
					} catch(NullPointerException e1) {
					}
				}
				if(numberTable.getSelectedRow()==2) {
					try {
					value = JOptionPane.showInputDialog(jf,"Enter new memory");
					if(value.length()>0) memory=Integer.valueOf(value);
					} catch(NullPointerException e1) {
					}
				}
				numberData[0][1]=numRanks;
				numberData[1][1]=numThreads;
				numberData[2][1]=memory;
				update();
			}
		});

		threshPanel = new JPanel();
		threshPanel.setLayout(new BoxLayout(threshPanel,BoxLayout.X_AXIS));
		threshPanel.setPreferredSize(new Dimension(140,70));
		threshPanel.setMinimumSize(new Dimension(140,70));
		threshPanel.setMaximumSize(new Dimension(140,70));
		LineBorder threshBorder = new LineBorder(Color.black);
		threshPanel.setBorder(threshBorder);
		Object[][] threshData = new Object[][] {
			{"thr_MO", thresh_MO},
			{"thr_CI", thresh_CI},
			{"IPEA", ipea}
		};
		String[] threshColumns = new String[] {" "," "};
		threshTable = new JTable(threshData,threshColumns);
		threshTable.setCellSelectionEnabled(true);
		threshTable.getColumnModel().getColumn(0).setMaxWidth(80);
		threshTable.getColumnModel().getColumn(1).setMaxWidth(50);
		ListSelectionModel threshSelectionModel = threshTable.getSelectionModel();
		threshSelectionModel.addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent e) {
				JFrame jf = new JFrame();
				String value;
				if(threshTable.getSelectedRow()==0) {
					try {
					value = JOptionPane.showInputDialog(jf,"Enter new threshold for common MO basis");
					if(value.length()>0) thresh_MO=Double.valueOf(value);
					} catch(NullPointerException e1) {
					}
				}
				if(threshTable.getSelectedRow()==1) {
					try {
					value = JOptionPane.showInputDialog(jf,"Enter new threshold for CI coefficients");
					if(value.length()>0) thresh_CI=Double.valueOf(value);
					} catch(NullPointerException e1) {
					}
				}
				if(threshTable.getSelectedRow()==2) {
					try {
					value = JOptionPane.showInputDialog(jf,"Enter new CASPT2 IPEA shift");
					if(value.length()>0) ipea=Double.valueOf(value);
					} catch(NullPointerException e1) {
					}
				}
				threshData[0][1]=thresh_MO;
				threshData[1][1]=thresh_CI;
				threshData[2][1]=ipea;
				update();
			}
		});
		
		expansionPanel = new JPanel();
		expansionPanel.setLayout(new BoxLayout(expansionPanel,BoxLayout.Y_AXIS));
		expansionPanel.setPreferredSize(new Dimension(140,70));
		expansionPanel.setMinimumSize(new Dimension(140,70));
		expansionPanel.setMaximumSize(new Dimension(140,70));
		LineBorder expansionBorder = new LineBorder(Color.black);
		expansionPanel.setBorder(expansionBorder);

		Object[][] expansionData = new Object[][] {
			{"Expansion", expMEBF[expansion]}
		};
		String[] expansionColumns = new String[] {" "," "};
		expansionTable = new JTable(expansionData,expansionColumns);
		expansionTable.setCellSelectionEnabled(true);
		expansionTable.getColumnModel().getColumn(0).setMaxWidth(80);
		expansionTable.getColumnModel().getColumn(1).setMaxWidth(100);

		
		DefaultComboBoxModel expansionComboModel = new DefaultComboBoxModel(expMEBF);
		JComboBox expansionCombo = new JComboBox();
		expansionCombo.setMaximumRowCount(3);
		expansionCombo.setModel(expansionComboModel);
		TableColumn expansionColumn = expansionTable.getColumnModel().getColumn(1);
		expansionColumn.setCellEditor(new DefaultCellEditor(expansionCombo));
		expansionCombo.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent item){
				Integer col = expansionTable.getEditingColumn();
				Integer row = expansionTable.getEditingRow();
				expansion=expansionCombo.getSelectedIndex();
//			    	                    if(row==1) contract=basisSetCombo.getSelectedIndex();
//			        	                if(row==2) cholesky=basisSetCombo.getSelectedIndex();
				update();
				for(int i=0; i<numMEBFs; i++) {
					Integer nmer = mebfSpecification[i][0];
					Integer spin = mebfSpecification[i][1];
					Integer chrg = mebfSpecification[i][2];
					Integer stat = mebfSpecification[i][3];
					selectMEBFStates(i,nmer,spin,chrg,stat);
				}
			}
		});

		
		basisPanel = new JPanel();
		basisPanel.setLayout(new BoxLayout(basisPanel,BoxLayout.Y_AXIS));
		basisPanel.setPreferredSize(new Dimension(140,70));
		basisPanel.setMinimumSize(new Dimension(140,70));
		basisPanel.setMaximumSize(new Dimension(140,70));
		LineBorder basisBorder = new LineBorder(Color.black);
		basisPanel.setBorder(basisBorder);
		
		Object[][] basisData = new Object[][] {
			{"Basis Set", basisSets[basisSet]}
		};
		String[] basisColumns = new String[] {" "," "};
		basisTable = new JTable(basisData,basisColumns);
		basisTable.setCellSelectionEnabled(true);
		basisTable.getColumnModel().getColumn(0).setMaxWidth(80);
		basisTable.getColumnModel().getColumn(1).setMaxWidth(100);
//		ListSelectionModel basisSelectionModel = basisTable.getSelectionModel();
		
//		ListSelectionModel basisSelectionModel = basisTable.getSelectionModel();
		
//		basisSelectionModel.addListSelectionListener(new ListSelectionListener() {
//			public void valueChanged(ListSelectionEvent e) {
				
		DefaultComboBoxModel basisSetComboModel = new DefaultComboBoxModel(basisSets);
		JComboBox basisSetCombo = new JComboBox();
		basisSetCombo.setMaximumRowCount(3);
		basisSetCombo.setModel(basisSetComboModel);
		TableColumn basisSetColumn = basisTable.getColumnModel().getColumn(1);
		basisSetColumn.setCellEditor(new DefaultCellEditor(basisSetCombo));
		basisSetCombo.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent item){
				Integer col = basisTable.getEditingColumn();
				Integer row = basisTable.getEditingRow();
				basisSet=basisSetCombo.getSelectedIndex();
//			    	                    if(row==1) contract=basisSetCombo.getSelectedIndex();
//			        	                if(row==2) cholesky=basisSetCombo.getSelectedIndex();
				update();
			}
		});

		/*
		contractPanel = new JPanel();
		contractPanel.setLayout(new BoxLayout(contractPanel,BoxLayout.X_AXIS));
		contractPanel.setPreferredSize(new Dimension(140,70));
		contractPanel.setMinimumSize(new Dimension(140,70));
		contractPanel.setMaximumSize(new Dimension(140,70));
		LineBorder contractBorder = new LineBorder(Color.black);
		contractPanel.setBorder(contractBorder);
		*/
		
		Object[][] contractData = new Object[][] {
			{"Contraction", contracts[contract]}
		};
		String[] contractColumns = new String[] {" "," "};
		contractTable = new JTable(contractData,contractColumns);
		contractTable.setCellSelectionEnabled(true);
		contractTable.getColumnModel().getColumn(0).setMaxWidth(80);
		contractTable.getColumnModel().getColumn(1).setMaxWidth(100);
		
		DefaultComboBoxModel contractComboModel = new DefaultComboBoxModel(contracts);
		JComboBox contractCombo = new JComboBox();
		contractCombo.setMaximumRowCount(3);
		contractCombo.setModel(contractComboModel);
		TableColumn contractColumn = contractTable.getColumnModel().getColumn(1);
		contractColumn.setCellEditor(new DefaultCellEditor(contractCombo));
		contractCombo.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent item){
				Integer col = contractTable.getEditingColumn();
				Integer row = contractTable.getEditingRow();
//			    	                    if(row==0) basisSet=basisSetCombo.getSelectedIndex();
				contract=contractCombo.getSelectedIndex();
//			    	                    if(row==2) cholesky=basisSetCombo.getSelectedIndex();
				update();
			}
		});	
				
/*
		choleskyPanel = new JPanel();
		choleskyPanel.setLayout(new BoxLayout(choleskyPanel,BoxLayout.X_AXIS));
		choleskyPanel.setPreferredSize(new Dimension(140,70));
		choleskyPanel.setMinimumSize(new Dimension(140,70));
		choleskyPanel.setMaximumSize(new Dimension(140,70));
		LineBorder choleskyBorder = new LineBorder(Color.black);
		choleskyPanel.setBorder(choleskyBorder);
		*/
		
		Object[][] choleskyData = new Object[][] {
			{"Cholesky", choleskys[cholesky]}
		};
		String[] choleskyColumns = new String[] {" "," "};
		choleskyTable = new JTable(choleskyData,choleskyColumns);
		choleskyTable.setCellSelectionEnabled(true);
		choleskyTable.getColumnModel().getColumn(0).setMaxWidth(80);
		choleskyTable.getColumnModel().getColumn(1).setMaxWidth(100);
				
		DefaultComboBoxModel choleskyComboModel = new DefaultComboBoxModel(choleskys);
		JComboBox choleskyCombo = new JComboBox();
		choleskyCombo.setMaximumRowCount(3);
		choleskyCombo.setModel(choleskyComboModel);
		TableColumn choleskyColumn = choleskyTable.getColumnModel().getColumn(1);
		choleskyColumn.setCellEditor(new DefaultCellEditor(choleskyCombo));
		choleskyCombo.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent item){
				Integer col = choleskyTable.getEditingColumn();
				Integer row = choleskyTable.getEditingRow();
//			    	                    if(row==0) basisSet=basisSetCombo.getSelectedIndex();
//			        	                if(row==1) contract=basisSetCombo.getSelectedIndex();
				cholesky=choleskyCombo.getSelectedIndex();
				update();
			}
		});	
		
		jobPanel = new JPanel();
		jobPanel.setLayout(new BoxLayout(jobPanel,BoxLayout.Y_AXIS));
		jobPanel.setPreferredSize(new Dimension(140,70));
		jobPanel.setMinimumSize(new Dimension(140,70));
		jobPanel.setMaximumSize(new Dimension(140,70));
		LineBorder jobBorder = new LineBorder(Color.black);
		jobPanel.setBorder(jobBorder);
		Object[][] jobData = new Object[][] {
			{"Account", account},
			{"Job Name", jobName},
			{"Time Limit", timeLimit}
		};
		String[] jobColumns = new String[] {" "," "};
		jobTable = new JTable(jobData,jobColumns);
		jobTable.setCellSelectionEnabled(true);
		jobTable.getColumnModel().getColumn(0).setMaxWidth(80);
		jobTable.getColumnModel().getColumn(1).setMaxWidth(100);
		jobTable.addMouseListener(new MouseAdapter() {
			public void mouseClicked(java.awt.event.MouseEvent e) {
				Integer row=jobTable.getSelectedRow();
				Integer col=jobTable.getSelectedColumn();
				JFrame jf = new JFrame();
				String value;
				if(row==0 && col==1) {
					try {
						value = JOptionPane.showInputDialog(jf,"Enter account");
						if(value.length()>0) account=value.trim();
						jobData[0][1]=account;
					} catch(NullPointerException e1) {
					}
				}
				if(row==1 && col==1) {
					try {
						value = JOptionPane.showInputDialog(jf,"Job name");
						if(value.length()>0) jobName=value.trim();
						jobData[1][1]=jobName;
					} catch(NullPointerException e1) {
					}
				}
				if(row==2 && col==1) {
					try {
						value = JOptionPane.showInputDialog(jf,"Enter time limit");
						if(value.length()>0) timeLimit=value.trim();
						jobData[2][1]=timeLimit;
					} catch(NullPointerException e1) {
					}
				}
				update();
			}
		});
		/*
        DefaultComboBoxModel basisSetComboModel = new DefaultComboBoxModel(basisSets);
        JComboBox basisSetCombo = new JComboBox();
        basisSetCombo.setMaximumRowCount(3);
        basisSetCombo.setModel(basisSetComboModel);
        
        TableCell basisSetCell = basisTable.getTableView();
//        TableColumn basisSetColumn = basisTable.getColumnModel().getColumn(1);
//        basisSetColumn.setCellEditor(new DefaultCellEditor(basisSetCombo));
        
        basisSetCombo.addItemListener(new ItemListener() {
            public void itemStateChanged(ItemEvent item){
                        Integer row = basisTable.getSelectedRow();
                        Integer col = basisTable.getSelectedColumn();
                        if(row==0) basisSet=basisSetCombo.getSelectedIndex();
                        if(row==1) contract=basisSetCombo.getSelectedIndex();
                        if(row==2) cholesky=basisSetCombo.getSelectedIndex();
                        update();
            }
        });
*/

		/*
		basisSelectionModel.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent item) {
				JFrame jf = new JFrame();
				String value;
				if(basisTable.getSelectedRow()==0) {
					try {
					value = JOptionPane.showInputDialog(jf,"Enter new threshold for common MO basis");
					if(value.length()>0) thresh_MO=Double.valueOf(value);
					} catch(NullPointerException e1) {
					}
				}
				if(basisTable.getSelectedRow()==1) {
					try {
					value = JOptionPane.showInputDialog(jf,"Enter new threshold for CI coefficients");
					if(value.length()>0) thresh_MO=Double.valueOf(value);
					} catch(NullPointerException e1) {
					}
				}
				basisData[0][1]=basisSets[basisSet];
				basisData[1][1]=contracts[contract];
				basisData[2][1]=choleskys[cholesky];
				update();
			}
		});
		*/
		
		JPanel buttonPanel = new JPanel();
		buttonPanel.setLayout(new BoxLayout(buttonPanel,BoxLayout.Y_AXIS));
		buttonPanel.setPreferredSize(new Dimension(115,70));
		buttonPanel.setMinimumSize(new Dimension(115,70));
		buttonPanel.setMaximumSize(new Dimension(115,70));
		LineBorder buttonBorder = new LineBorder(Color.black);
		buttonPanel.setBorder(buttonBorder);
		JButton updateButton = new JButton("Update");
		updateButton.setPreferredSize(new Dimension(110,20));
		updateButton.setMinimumSize(new Dimension(110,20));
		updateButton.setMaximumSize(new Dimension(110,20));
		updateButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				update();
				System.out.println("Data updated");
			}
		});

		JButton generateButton = new JButton("Generate");
		generateButton.setPreferredSize(new Dimension(110,20));
		generateButton.setMinimumSize(new Dimension(110,20));
		generateButton.setMaximumSize(new Dimension(110,20));
		generateButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				update();
				writeInputFiles();
				writeClearScripts();
				System.out.println("Files generated");
			}
		});
		
		optionsPanel.add(Box.createRigidArea(new Dimension(5,5)));
		optionsPanel.add(dimersButton);
		optionsPanel.add(dimersNoAButton);
		optionsPanel.add(Box.createVerticalGlue());
		dimensionPanel.add(Box.createRigidArea(new Dimension(5,0)));
		dimensionPanel.add(dimensionTable);
		numberPanel.add(Box.createRigidArea(new Dimension(5,0)));
		numberPanel.add(numberTable);
		threshPanel.add(Box.createRigidArea(new Dimension(5,0)));
		threshPanel.add(threshTable);
		basisPanel.add(Box.createRigidArea(new Dimension(5,0)));
		basisPanel.add(basisTable);
		basisPanel.add(contractTable);
		basisPanel.add(choleskyTable);
		expansionPanel.add(Box.createRigidArea(new Dimension(5,0)));
		expansionPanel.add(expansionTable);
		jobPanel.add(Box.createRigidArea(new Dimension(5,0)));
		jobPanel.add(jobTable);
		buttonPanel.add(Box.createRigidArea(new Dimension(5,5)));
		buttonPanel.add(updateButton);
		buttonPanel.add(generateButton);
		buttonPanel.add(Box.createVerticalGlue());
		parametersPanel.add(dimensionPanel);
		parametersPanel.add(Box.createRigidArea(new Dimension(5,0)));
		parametersPanel.add(optionsPanel);
		parametersPanel.add(Box.createRigidArea(new Dimension(5,0)));
		parametersPanel.add(numberPanel);
		parametersPanel.add(Box.createRigidArea(new Dimension(5,0)));
		parametersPanel.add(threshPanel);
		parametersPanel.add(Box.createRigidArea(new Dimension(5,0)));
		parametersPanel.add(expansionPanel);
		parametersPanel.add(Box.createRigidArea(new Dimension(5,0)));
		parametersPanel.add(basisPanel);
		parametersPanel.add(Box.createRigidArea(new Dimension(5,0)));
		parametersPanel.add(jobPanel);
		parametersPanel.add(Box.createRigidArea(new Dimension(5,0)));
		parametersPanel.add(buttonPanel);
		
		parametersPanel.add(Box.createHorizontalGlue());
		baseBox.add(parametersPanel);
		
// States Panel

		statesPanel = new JPanel();
		statesPanel.setLayout(new BoxLayout(statesPanel,BoxLayout.Y_AXIS));
		statesPanel.setPreferredSize(new Dimension(Short.MAX_VALUE,((numSets)*15+35)));
		statesPanel.setMinimumSize(new Dimension(Short.MAX_VALUE,((numSets)*15+35)));
		TitledBorder statesBorder = new TitledBorder(new LineBorder(Color.black),"State List Definitions");
		statesBorder.setTitleColor(Color.black);
		statesPanel.setBorder(statesBorder);
		statesTable = new JTable(statesTableModel);
		statesTable.getColumnModel().getColumn(0).setMaxWidth(10);
		statesTable.getColumnModel().getColumn(1).setMaxWidth(40);
		statesTable.getColumnModel().getColumn(2).setMaxWidth(40);
		ListSelectionModel statesSelectionModel = statesTable.getSelectionModel();
		statesSelectionModel.addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent e) {
				Integer row = statesTable.getSelectedRow();
				Integer col = statesTable.getSelectedColumn();
				update();
			}
		});
		JScrollPane statesScroll = new JScrollPane(statesTable);
		statesTable.setCellSelectionEnabled(true);
		statesTable.addMouseListener(new MouseAdapter() {
			public void mouseReleased(MouseEvent e) {
					stateCellSelected(e);
			}
		});
		
		DefaultComboBoxModel setComboModel = new DefaultComboBoxModel(stateNames);
		JComboBox setCombo = new JComboBox();
		setCombo.setMaximumRowCount(numStateNames);
		setCombo.setModel(setComboModel);
		for(int i=0; i<14; i++) {
			TableColumn setColumn = statesTable.getColumnModel().getColumn(3+i);
			setColumn.setCellEditor(new DefaultCellEditor(setCombo));
		}
		setCombo.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent item) {
				Integer row = statesTable.getEditingRow();
				Integer col = statesTable.getEditingColumn();
				if(row>=0 && col>=3 && item.getStateChange()==ItemEvent.SELECTED) {
					ndxStateList[row][(col-3)]=setCombo.getSelectedIndex();
					update();
				}
			}
		});
		
		statesPanel.add(statesScroll);
		baseBox.add(statesPanel);

// Fragment Panel
		
		fragmentsPanel = new JPanel();
		fragmentsPanel.setLayout(new BoxLayout(fragmentsPanel,BoxLayout.Y_AXIS));
		fragmentsPanel.setPreferredSize(new Dimension(Short.MAX_VALUE,((numFragments)*15+35)));
		fragmentsPanel.setMinimumSize(new Dimension(Short.MAX_VALUE,((numFragments)*15+35)));
		TitledBorder fragmentsBorder = new TitledBorder(new LineBorder(Color.black),"Fragment Definitions");
		fragmentsBorder.setTitleColor(Color.black);
		fragmentsPanel.setBorder(fragmentsBorder);
		fragmentsTable = new JTable(fragmentsTableModel);
		fragmentsTable.getColumnModel().getColumn(0).setMaxWidth(10);
		fragmentsTable.getColumnModel().getColumn(1).setMaxWidth(30);
		fragmentsTable.getColumnModel().getColumn(3).setMaxWidth(50);
		fragmentsTable.getColumnModel().getColumn(4).setMaxWidth(50);
		fragmentsTable.getColumnModel().getColumn(5).setMaxWidth(70);
		fragmentsTable.getColumnModel().getColumn(6).setMaxWidth(50);
		fragmentsTable.getColumnModel().getColumn(7).setMaxWidth(50);
		fragmentsTable.getColumnModel().getColumn(8).setMaxWidth(50);
		fragmentsTable.getColumnModel().getColumn(9).setMaxWidth(50);
		fragmentsTable.getColumnModel().getColumn(10).setMaxWidth(50);
		fragmentsTable.getColumnModel().getColumn(11).setMaxWidth(50);
		fragmentsTable.getColumnModel().getColumn(12).setMaxWidth(50);
		fragmentsTable.getColumnModel().getColumn(13).setMaxWidth(50);
		fragmentsTable.getColumnModel().getColumn(14).setMaxWidth(50);
		
		
//		ListSelectionModel fragmentSelectionModel = fragmentsTable.getSelectionModel();
//		fragmentSelectionModel.addListSelectionListener(new ListSelectionListener() {
//			public void valueChanged(ListSelectionEvent e) {

		fragmentsTable.addMouseListener(new MouseAdapter() {
			public void mouseClicked(java.awt.event.MouseEvent e) {			
				
				Integer row = fragmentsTable.getSelectedRow();
				Integer col = fragmentsTable.getSelectedColumn();
				JFrame jf = new JFrame();
				String value;
				String nameA, nameB, nameP, nameF;
				
				// xyz coordinate file
				if(col==2) {
					for(int j=0; j<numFragments; j++) {
					}
					String newFile = getFileRoot(".xyz");
					if(newFile.length()>0) {
						fragmentDefinitions[row][col]=newFile.trim();
					}
					namFragments[row]=fragmentDefinitions[row][col].toString();
					dimFragments[row][1]=getNumAxyz(fragmentDefinitions[row][col].toString());
					dimFragments[row][3]=getNumExyz(fragmentDefinitions[row][col].toString());
					Integer numXA=dimFragments[row][1]-getNumHxyz(fragmentDefinitions[row][col].toString());
					if(numXA<dimFragments[row][4]) {
						numCASe=Math.max(4,numXA);
						numCASo=Math.max(4,numXA);
						dimFragments[row][4]=numCASe;
						dimFragments[row][5]=numCASo;
					}
					for(int j=0; j<numFragments; j++) {
						if(j!=row) {
							if(namFragments[j].length()>0) {
								dimFragments[row][0]=row;
								for(int i=0; i<numFragments; i++) {
									if(i!=row && namFragments[i]==namFragments[row]) {
										if(i<row) {
											dimFragments[row][0]=i;
										} else {
											dimFragments[i][0]=row;
										}
									}
								} 
							} else {
								if(fragmentDefinitions[j][1]==fragmentDefinitions[row][1]) {
									dimFragments[j][1]=dimFragments[row][1];
									dimFragments[j][3]=dimFragments[row][3];
									dimFragments[j][4]=dimFragments[row][4];
									dimFragments[j][5]=dimFragments[row][5];
									namFragments[j]=namFragments[row];
								}
							}
						}
					}
					update();
				}
				// change number of states from { S0 S1 T1 D- D+ S2 T2 }
				if(col==4) {
					try {
						value = JOptionPane.showInputDialog(jf,"Enter index into states list for fragment "+fragmentDefinitions[row][0].toString());
						if(value.length()>0) dimFragments[row][2]=Integer.valueOf(value)-1;
					} catch(NullPointerException e1) {
					}
				}
				// change number of CAS electrons
				if(col==6) {
					try {
						value = JOptionPane.showInputDialog(jf,"Enter number of electrons in CAS for fragment "+fragmentDefinitions[row][0].toString());
						if(value.length()>0) dimFragments[row][4]=Integer.valueOf(value);
					} catch(NullPointerException e1) {
					}
				}
				// change number of CAS orbitals
				if(col==7) {
					try {
						value = JOptionPane.showInputDialog(jf,"Enter number of orbitals in CAS for fragment "+fragmentDefinitions[row][0].toString());
						if(value.length()>0) dimFragments[row][5]=Integer.valueOf(value);
					} catch(NullPointerException e1) {
					}
				}
				// change translation in x
				if(col==8) {
					try {
						value = JOptionPane.showInputDialog(jf,"Enter Tx for fragment "+fragmentDefinitions[row][0].toString());
						movFragments[row][0]=Double.valueOf(value).doubleValue();
					} catch(NullPointerException e1) {
					}
				}
				// change translation in y
				if(col==9) {
					try {
						value = JOptionPane.showInputDialog(jf,"Enter Ty for fragment "+fragmentDefinitions[row][0].toString());
						movFragments[row][1]=Double.valueOf(value).doubleValue();
					} catch(NullPointerException e1) {
					}
				}
				// change translation in z
				if(col==10) {
					try {
						value = JOptionPane.showInputDialog(jf,"Enter Tz for fragment "+fragmentDefinitions[row][0].toString());
						movFragments[row][2]=Double.valueOf(value).doubleValue();
					} catch(NullPointerException e1) {
					}
				}
				// change rotation in x
				if(col==11) {
					try {
						value = JOptionPane.showInputDialog(jf,"Enter Rx for fragment "+fragmentDefinitions[row][0].toString());
						movFragments[row][3]=Double.valueOf(value).doubleValue();
					} catch(NullPointerException e1) {
					}
				}
				// change rotation in y
				if(col==12) {
					try {
						value = JOptionPane.showInputDialog(jf,"Enter Ry for fragment "+fragmentDefinitions[row][0].toString());
						movFragments[row][4]=Double.valueOf(value).doubleValue();
					} catch(NullPointerException e1) {
					}
				}
				// change rotation in z
				if(col==13) {
					try {
						value = JOptionPane.showInputDialog(jf,"Enter Rz for fragment "+fragmentDefinitions[row][0].toString());
						movFragments[row][5]=Double.valueOf(value).doubleValue();
					} catch(NullPointerException e1) {
					}
				}
				if(col==14) {
					for(int j=0; j<6; j++) RandT[j]=movFragments[row][j];
					nameP=projectName.trim();
					nameF=namFragments[row].trim();
					nameA= (String) fragmentDefinitions[dimFragments[row][0]][0];
					nameB= (String) namFragments[row];
					fragment.initialize(nameP,nameF, nameA,nameB,dimFragments[row][3],RandT);
				}
				update();
			}
		});

		
		DefaultComboBoxModel equivalenceComboModel = new DefaultComboBoxModel(fragmentNames);
		JComboBox equivalenceCombo = new JComboBox();
		equivalenceCombo.setMaximumRowCount(numFragments);
		equivalenceCombo.setModel(equivalenceComboModel);
		TableColumn equivalenceColumn = fragmentsTable.getColumnModel().getColumn(1);
		equivalenceColumn.setCellEditor(new DefaultCellEditor(equivalenceCombo));
		equivalenceCombo.addItemListener(new ItemListener() {
		    public void itemStateChanged(ItemEvent item){
				Integer row = fragmentsTable.getEditingRow();
				Integer col = fragmentsTable.getEditingColumn();
				fragmentDefinitions[row][col]=fragmentNames[equivalenceCombo.getSelectedIndex()];
				update();
		    }
		});
		
		JScrollPane fragmentScroll = new JScrollPane(fragmentsTable);
		fragmentsPanel.add(fragmentScroll);
		baseBox.add(fragmentsPanel);
		
// Energies Panel
		
		energiesPanel = new JPanel();
		energiesPanel.setLayout(new BoxLayout(energiesPanel,BoxLayout.X_AXIS));
		TitledBorder energiesBorder = new TitledBorder(new LineBorder(Color.black),"Fragment Energies");
		energiesBorder.setTitleColor(Color.black);
		energiesPanel.setBorder(energiesBorder);
		energiesTable = new JTable(energiesTableModel);
		JScrollPane energiesScroll = new JScrollPane(energiesTable);
		energiesTable.setCellSelectionEnabled(true);
		energiesTable.addMouseListener(new MouseAdapter() {
			public void mouseReleased(MouseEvent e) {
					energyCellSelected(e);
			}
		});
		energiesPanel.add(energiesScroll);
		baseBox.add(energiesPanel);
				
		mebfsPanel = new JPanel();
		mebfsPanel.setLayout(new BoxLayout(mebfsPanel,BoxLayout.X_AXIS));
		mebfsPanel.setPreferredSize(new Dimension(Short.MAX_VALUE,numberMebfDefinitions*15+35));
		mebfsPanel.setMinimumSize(new Dimension(Short.MAX_VALUE,numberMebfDefinitions*15+35));
		TitledBorder mebfsBorder = new TitledBorder(new LineBorder(Color.black),"MEBF List Definitions");
		mebfsBorder.setTitleColor(Color.black);
		mebfsPanel.setBorder(mebfsBorder);
		mebfsTable = new JTable(mebfsTableModel);
		mebfsTable.getColumnModel().getColumn(0).setMaxWidth(60);
		for(int i=6; i<25; i++) mebfsTable.getColumnModel().getColumn(i).setMaxWidth(60);
		ListSelectionModel mebfsSelectionModel = mebfsTable.getSelectionModel();
		mebfsSelectionModel.addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent e) {
				Integer row = mebfsTable.getSelectedRow();
				Integer col = mebfsTable.getSelectedColumn();
				update();
			}
		});
		JScrollPane mebfsScroll = new JScrollPane(mebfsTable);
		mebfsTable.setCellSelectionEnabled(true);
		mebfsTable.addMouseListener(new MouseAdapter() {
			public void mouseReleased(MouseEvent e) {
					mebfCellSelected(e);
			}
		});

		DefaultComboBoxModel fragmentComboModel = new DefaultComboBoxModel(fragmentNames);
		JComboBox fragmentCombo = new JComboBox();
		fragmentCombo.setModel(fragmentComboModel);
		TableColumn fragmentColumn = mebfsTable.getColumnModel().getColumn(5);
		fragmentColumn.setCellEditor(new DefaultCellEditor(fragmentCombo));
		fragmentCombo.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent item) {
				Integer row = mebfsTable.getEditingRow();
				Integer col = mebfsTable.getEditingColumn();
				if(row>=0 && col==5) {
					Integer mymebf = mebfIndex[row][0];
					Integer mynmer = mebfIndex[row][1];
					mebfFragments[mymebf][mynmer][0] = fragmentCombo.getSelectedIndex();
				}
			}
		});
		
		DefaultComboBoxModel stateComboModel = new DefaultComboBoxModel(stateNames);
		JComboBox stateCombo = new JComboBox();
		stateCombo.setModel(stateComboModel);
		for(int i=0; i<maxMEBFstates; i++) {
			TableColumn stateColumn = mebfsTable.getColumnModel().getColumn(6+i);
			stateColumn.setCellEditor(new DefaultCellEditor(stateCombo));
		}
		stateCombo.addItemListener(new ItemListener() {
			public void itemStateChanged(ItemEvent item) {
				Integer row = mebfsTable.getEditingRow();
				Integer col = mebfsTable.getEditingColumn();
				if(row>=0 && col>=5) {
					Integer mymebf = mebfIndex[row][0];
					Integer mynmer = mebfIndex[row][1];
					mebfFragments[mymebf][mynmer][col-4] = stateCombo.getSelectedIndex();
					mebfDefinitions[row][col]=stateNames[mebfFragments[mymebf][mynmer][col-4]];
					mebfsPanel.revalidate();
				}
				update();
			}
		});
		
		mebfsPanel.add(mebfsScroll);
		baseBox.add(mebfsPanel);

		nociResultsPanel = new JPanel();
		nociResultsPanel.setLayout(new BoxLayout(nociResultsPanel,BoxLayout.X_AXIS));
		TitledBorder nociResultsBorder = new TitledBorder(new LineBorder(Color.black),"NOCI Energies");
		nociResultsBorder.setTitleColor(Color.black);
		nociResultsPanel.setBorder(nociResultsBorder);
		nociResultsPanel.setPreferredSize(new Dimension(Short.MAX_VALUE,noci1*20+135));
		nociResultsPanel.setMinimumSize(new Dimension(Short.MAX_VALUE,noci1*20+135));
		nociResultsPanel.setMaximumSize(new Dimension(Short.MAX_VALUE,noci1*20+135));
		
		nociEnergiesPanel = new JPanel();
		nociEnergiesPanel.setLayout(new BoxLayout(nociEnergiesPanel,BoxLayout.X_AXIS));
		nociEnergiesPanel.setPreferredSize(new Dimension(noci0*80+40,noci1*20+130));
		nociEnergiesPanel.setMinimumSize(new Dimension(noci0*80+40,noci1*20+80));
		nociEnergiesPanel.setMaximumSize(new Dimension(noci0*80+40,noci1*20+80));
		LineBorder nociEnergiesBorder = new LineBorder(Color.black);
		nociEnergiesPanel.setBorder(nociEnergiesBorder);
		nociEnergiesTableModel = new DefaultTableModel(new Object[] {"ID"},1);
		noci0=1;
		noci1=1;
		nociEnergiesTable = new JTable(nociEnergiesTableModel);
		nociEnergiesTable.getColumnModel().getColumn(0).setMaxWidth(40);
		ListSelectionModel nociEnergiesSelectionModel = nociEnergiesTable.getSelectionModel();
		nociEnergiesSelectionModel.addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent e) {
				Integer row = nociEnergiesTable.getSelectedRow();
				Integer col = nociEnergiesTable.getSelectedColumn();
				nociEnergiesPanel.revalidate();
				nociEnergiesPanel.repaint();
			}
		});
		JScrollPane nociEnergiesScroll = new JScrollPane(nociEnergiesTable);
		nociEnergiesTable.setCellSelectionEnabled(true);
		nociEnergiesTable.addMouseListener(new MouseAdapter() {
			public void mouseReleased(MouseEvent e) {
					nociEnergiesCellSelected(e);
			}
		});

		nociEnergiesPanel.add(nociEnergiesScroll);

		nociPlotPanel = new JPanel();
		nociPlotPanel.setLayout(new BoxLayout(nociPlotPanel,BoxLayout.X_AXIS));
		nociPlotPanel.setPreferredSize(new Dimension(300,noci1*20+130));
		nociPlotPanel.setMinimumSize(new Dimension(300,800));
		nociPlotPanel.setMaximumSize(new Dimension(300,800));
		LineBorder nociPlotBorder = new LineBorder(Color.black);
		nociPlotPanel.setBorder(nociPlotBorder);
		energyDiagram.setXRange(0,3);
		energyDiagram.setYRange(0,10);
		energyDiagram.setSize(200,500);
		energyDiagram.setGrid(false);
		nociPlotPanel.add(energyDiagram);
		
		nociResultsPanel.add(nociEnergiesPanel);
		nociResultsPanel.add(Box.createRigidArea(new Dimension(5,0)));
		nociResultsPanel.add(nociPlotPanel);
		nociResultsPanel.add(Box.createHorizontalGlue());
		nociResultsPanel.add(Box.createVerticalGlue());

		energyDiagram.revalidate();
		energyDiagram.repaint();
		nociPlotPanel.revalidate();
		nociPlotPanel.repaint();
		
		baseBox.add(nociResultsPanel);
		baseBox.add(Box.createVerticalGlue());

		update();
	}

	private void stateCellSelected(MouseEvent e) {
		Integer row = statesTable.getSelectedRow();
		Integer col = statesTable.getSelectedColumn();
		JFrame jf = new JFrame();
		String value;
		if(col==1) {
			try {
				Integer prev = lenStateList[row];
				value = JOptionPane.showInputDialog(jf,"Enter number of states for list "+row);
				if(value.length()>0) lenStateList[row]=Integer.valueOf(value.trim());
				if(prev<lenStateList[row]) {
					for(int i=prev; i<lenStateList[row]; i++) ndxStateList[row][i]=0;
				}
			} catch (NullPointerException e1) {
			}
		}
		if(col==2) {
			try {
				value = JOptionPane.showInputDialog(jf,"Enter spin of states for list "+row);
				if(value.length()>0) spinStateList[row]=Integer.valueOf(value.trim());
			} catch (NullPointerException e1) {
			}
		}
		updateStatesList();
	}
	
	private void energyCellSelected(MouseEvent e) {
	}


	private void selectMEBFStates(Integer mebf, Integer nmer, Integer spin, Integer charge, Integer nums) {

		Integer prev=nums;
		Integer curr=mebfSpecification[mebf][3];
		Integer count = 0;
		Integer[] lens = new Integer[5];
		Integer[] ndxs = new Integer[5];
		String[] sw1 = new String[] {"S","D","T","q","Q"};
		String[] sw2 = new String[] {"D","S","T","q","Q"};
		String[][] sw = new String[5][5];
		Boolean include = false;
		
		if(mebfSpecification[mebf][0]>maxMer) mebfSpecification[mebf][0]=maxMer;
		
		if(mebfSpecification[mebf][0]>nmer) {
			for(int i=nmer; i<mebfSpecification[mebf][0]; i++) {
				for(int k=0; k<maxMEBFstates; k++) {
					mebfFragments[mebf][i][k]=0;
				}
				for(int k=0; k<mebfSpecification[mebf][4]+1; k++) {
					mebfFragments[mebf][i][k]=mebfFragments[mebf][nmer-1][k];
				}
				mebfFragments[mebf][i][0]=mebfFragments[mebf][i-1][0]+1;
				if(mebfFragments[mebf][i][0]>maxMer) mebfFragments[mebf][i][0]=mebfFragments[mebf][i][0]-maxMer;
				
			}	
		}
		
		if(mebfSpecification[mebf][0]==1 && numSets==1) mebfSpecification[mebf][1]=spinStateList[0];
					
		nmer=mebfSpecification[mebf][0];
		
		spin=mebfSpecification[mebf][1];

		charge=mebfSpecification[mebf][2];
		
// Initialize monomers
//		if(numMEBFs==3 && numFragments==3 &&
//					mebfSpecification[0][0]==1 && mebfSpecification[mebf][0]==1 && 
//					mebfSpecification[0][1]==1 && mebfSpecification[1][1]==1 && mebfSpecification[2][1]==1 && 
//					mebfSpecification[0][4]==1 && mebfSpecification[1][4]==1 && mebfSpecification[2][4]==1) {
//			if(mebf==1 && mebfFragments[0][0][0]==0 && mebfFragments[1][0][0]==0) mebfFragments[1][0][0]=1;
//			if(mebf==2 && mebfFragments[0][0][0]==0 && mebfFragments[2][0][0]==0) mebfFragments[2][0][0]=2;
//		}
		
//      Initialize dimers		
		if(numMEBFs==3 && numFragments==3 &&
					mebfSpecification[0][0]==2 && mebfSpecification[mebf][0]==2 && 
					mebfSpecification[0][1]==1 && mebfSpecification[1][1]==1 && mebfSpecification[2][1]==1 && 
					mebfSpecification[0][4]==1 && mebfSpecification[1][4]==1 && mebfSpecification[2][4]==1) {
			if(mebf==1 && mebfFragments[0][0][0]==0 && mebfFragments[0][1][0]==1 && 
					mebfFragments[1][0][0]==0 && mebfFragments[1][1][0]==1 ) {
				mebfFragments[1][0][0]=1;
				mebfFragments[1][1][0]=2;
			}
			if(mebf==2 && mebfFragments[0][0][0]==0 && mebfFragments[0][1][0]==1 && 
					mebfFragments[2][0][0]==0 && mebfFragments[2][1][0]==1 ) {
				mebfFragments[2][0][0]=0;
				mebfFragments[2][1][0]=2;
			}
		}
		
		if(numMEBFs==3 && numFragments==2 &&
				mebfSpecification[0][0]==2 && mebfSpecification[mebf][0]==2 && 
				mebfSpecification[0][1]==1 && mebfSpecification[1][1]==1 && mebfSpecification[2][1]==1 && 
				mebfSpecification[0][4]==1 && mebfSpecification[1][4]==1 && mebfSpecification[2][4]==1) {
			if(mebf==1 && mebfFragments[0][0][0]==0 && mebfFragments[0][1][0]==1 && 
					mebfFragments[1][0][0]==0 && mebfFragments[1][1][0]==1 ) {
				mebfSpecification[mebf][1]=3;
			}
			if(mebf==2 && mebfFragments[0][0][0]==0 && mebfFragments[0][1][0]==1 && 
					mebfFragments[2][0][0]==0 && mebfFragments[2][1][0]==1 ) {
				mebfSpecification[mebf][1]=5;
			}
		}

		// For monomer with spin 1
		if(mebfSpecification[mebf][0]==1 && mebfSpecification[mebf][1]==1) {
			Integer stateSet = dimFragments[mebfFragments[mebf][0][0]][2];
			Integer stateLen = lenStateList[stateSet];
			if(curr<prev) {
				mebfSpecification[mebf][3]=curr;
			} else {
				count=0;
				for(int i=0; i<stateLen; i++) {
					include=stateNames[ndxStateList[stateSet][i]].trim().startsWith("S");
					if(charge==0 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("+")>=0) include=false;
					if(charge==0 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("-")>=0) include=false;
					if(charge==1 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("+")<0) include=false;
					if(charge==-1 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("-")<0) include=false;
					if(include) {
						mebfFragments[mebf][0][count+1]=ndxStateList[stateSet][i];
						count++;
					}
				}
				mebfSpecification[mebf][3]=count;
			}
		}
		// For monomer with spin 2
		if(mebfSpecification[mebf][0]==1 && mebfSpecification[mebf][1]==2) {
			Integer stateSet = dimFragments[mebfFragments[mebf][0][0]][2];
			Integer stateLen = lenStateList[stateSet];
			if(curr<prev) {
				mebfSpecification[mebf][3]=curr;
			} else {
				count=0;
				for(int i=0; i<stateLen; i++) {
					include=false;
					if(stateNames[ndxStateList[stateSet][i]].trim().startsWith("D")) {
						include=true;
						if(charge==0 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("+")>=0) include=false;
						if(charge==0 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("-")>=0) include=false;
						if(charge==1 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("+")<0) include=false;
						if(charge==-1 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("-")<0) include=false;
					}
					if(stateNames[ndxStateList[stateSet][i]].trim().startsWith("S")) {
						include=false;
						if(charge==0 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("+")>=0) include=true;
						if(charge==0 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("-")>=0) include=true;
						if(charge==1 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("+")<0) include=true;
						if(charge==-1 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("-")<0) include=true;
					}
					if(stateNames[ndxStateList[stateSet][i]].trim().startsWith("T")) {
						include=false;
						if(charge==0 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("+")>=0) include=true;
						if(charge==0 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("-")>=0) include=true;
						if(charge==1 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("+")<0) include=true;
						if(charge==-1 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("-")<0) include=true;
					}
					if(include) {
						mebfFragments[mebf][0][count+1]=ndxStateList[stateSet][i];
						count++;
					}
				}
				mebfSpecification[mebf][3]=count;
			}

		}
			// For monomer with spin 3
		if(mebfSpecification[mebf][0]==1 && mebfSpecification[mebf][1]==3) {
			Integer stateSet = dimFragments[mebfFragments[mebf][0][0]][2];
			Integer stateLen = lenStateList[stateSet];
			if(curr<prev) {
				mebfSpecification[mebf][3]=curr;
			} else {
				count=0;
				for(int i=0; i<stateLen; i++) {
					include=stateNames[ndxStateList[stateSet][i]].trim().startsWith("T");
					if(charge==0 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("+")>=0) include=false;
					if(charge==0 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("-")>=0) include=false;
					if(charge==1 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("+")<0) include=false;
					if(charge==-1 && stateNames[ndxStateList[stateSet][i]].trim().indexOf("-")<0) include=false;
					if(include) {
						mebfFragments[mebf][0][count+1]=ndxStateList[stateSet][i];
						count++;
					}
				}
				mebfSpecification[mebf][3]=count;
			}
		}

				// For dimer with spin 1
		if(mebfSpecification[mebf][0]==2 && mebfSpecification[mebf][1]==1) {
			for(int k=0; k<2; k++) {
				ndxs[k]=dimFragments[mebfFragments[mebf][k][0]][2];
				lens[k]=lenStateList[ndxs[k]];
				if(spinStateList[ndxs[k]]==1 || spinStateList[ndxs[k]]==3 || spinStateList[ndxs[k]]==5) {
					for(int l=0; l<5; l++) sw[k][l]=sw1[l];
				} else {
					for(int l=0; l<5; l++) sw[k][l]=sw2[l];
				}
			}
			if(curr<prev) {
				mebfSpecification[mebf][3]=curr;
			} else {
				count=0;
				for(int i=0; i<5; i++) {
					for(int j0=0; j0<lens[0]; j0++) {
						for(int j1=0; j1<lens[0]; j1++) {
							include=stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith(sw[0][i]) && stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith(sw[1][i]);
							if(charge==0 && stateNames[ndxStateList[ndxs[0]][j0]].trim().indexOf("+")>=0 && stateNames[ndxStateList[ndxs[1]][j1]].trim().indexOf("-")<0) include=false;
							if(charge==0 && stateNames[ndxStateList[ndxs[0]][j0]].trim().indexOf("-")>=0 && stateNames[ndxStateList[ndxs[1]][j1]].trim().indexOf("+")<0) include=false;
							if(include) {
								mebfFragments[mebf][0][count+1]=ndxStateList[ndxs[0]][j0];
								mebfFragments[mebf][1][count+1]=ndxStateList[ndxs[1]][j1];
								count++;
							}
						}
					}
				}
				mebfSpecification[mebf][3]=count;
			}
		}

		// For dimer with spin 2
		if(mebfSpecification[mebf][0]==2 && mebfSpecification[mebf][1]==2) {
			for(int k=0; k<2; k++) {
				ndxs[k]=dimFragments[mebfFragments[mebf][k][0]][2];
				lens[k]=lenStateList[ndxs[k]];
				if(spinStateList[ndxs[k]]==1 || spinStateList[ndxs[k]]==3) {
					for(int l=0; l<3; l++) sw[k][l]=sw1[l];
				} else {
					for(int l=0; l<3; l++) sw[k][l]=sw2[l];
				}
			}
			if(curr<prev) {
				mebfSpecification[mebf][3]=curr;
			} else {
				count=0;
				for(int j0=0; j0<lens[0]; j0++) {
					for(int j1=0; j1<lens[0]; j1++) {					
						include=false;
						if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("S") && stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("D")) include=true;
						if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("D") && stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("S")) include=true;
						if(include) {
							mebfFragments[mebf][0][count+1]=ndxStateList[ndxs[0]][j0];
							mebfFragments[mebf][1][count+1]=ndxStateList[ndxs[1]][j1];
							count++;
						}
					}
				}
				mebfSpecification[mebf][3]=count;
			}
		}

		// For dimer with spin 3
		if(mebfSpecification[mebf][0]==2 && mebfSpecification[mebf][1]==3) {
			for(int k=0; k<2; k++) {
				ndxs[k]=dimFragments[mebfFragments[mebf][k][0]][2];
				lens[k]=lenStateList[ndxs[k]];
				if(spinStateList[ndxs[k]]==1 || spinStateList[ndxs[k]]==3) {
					for(int l=0; l<3; l++) sw[k][l]=sw1[l];
				} else {
					for(int l=0; l<3; l++) sw[k][l]=sw2[l];
				}
			}
			if(curr<prev) {
				mebfSpecification[mebf][3]=curr;
			} else {
				count=0;
				for(int j0=0; j0<lens[0]; j0++) {
					for(int j1=0; j1<lens[0]; j1++) {
						include=false;
						if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith(sw[0][0]) && stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith(sw[1][2])) include=true;
						if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith(sw[0][2]) && stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith(sw[1][0])) include=true;
						if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith(sw[0][2]) && stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith(sw[1][2])) include=true;
						if(charge==0 && stateNames[ndxStateList[ndxs[0]][j0]].trim().indexOf("+")>=0 && stateNames[ndxStateList[ndxs[1]][j1]].trim().indexOf("-")<0) include=false;
						if(charge==0 && stateNames[ndxStateList[ndxs[0]][j0]].trim().indexOf("-")>=0 && stateNames[ndxStateList[ndxs[1]][j1]].trim().indexOf("+")<0) include=false;
						if(include) {
							mebfFragments[mebf][0][count+1]=ndxStateList[ndxs[0]][j0];
							mebfFragments[mebf][1][count+1]=ndxStateList[ndxs[1]][j1];
							count++;
						}
					}
				}
				mebfSpecification[mebf][3]=count;
			}	
		}

		// For dimer with spin 5
		if(mebfSpecification[mebf][0]==2 && mebfSpecification[mebf][1]==5) {
			for(int k=0; k<2; k++) {
				ndxs[k]=dimFragments[mebfFragments[mebf][k][0]][2];
				lens[k]=lenStateList[ndxs[k]];
				if(spinStateList[ndxs[k]]==1 || spinStateList[ndxs[k]]==3) {
					for(int l=0; l<3; l++) sw[k][l]=sw1[l];
				} else {
					for(int l=0; l<3; l++) sw[k][l]=sw2[l];
				}
			}
			if(curr<prev) {
				mebfSpecification[mebf][3]=curr;
			} else {
				count=0;
				for(int j0=0; j0<lens[0]; j0++) {
					for(int j1=0; j1<lens[0]; j1++) {
						include=false;
						if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("T") && stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("T")) include=true;
						if(charge==0 && stateNames[ndxStateList[ndxs[0]][j0]].trim().indexOf("+")>=0 && stateNames[ndxStateList[ndxs[1]][j1]].trim().indexOf("-")<0) include=false;
						if(charge==0 && stateNames[ndxStateList[ndxs[0]][j0]].trim().indexOf("-")>=0 && stateNames[ndxStateList[ndxs[1]][j1]].trim().indexOf("+")<0) include=false;
						if(include) {
							mebfFragments[mebf][0][count+1]=ndxStateList[ndxs[0]][j0];
							mebfFragments[mebf][1][count+1]=ndxStateList[ndxs[1]][j1];
							count++;
						}
					}
				}
				mebfSpecification[mebf][3]=count;
			}	
		}
		
		// For trimer with spin 1
		if(mebfSpecification[mebf][0]==3 && mebfSpecification[mebf][1]==1) {
			for(int k=0; k<3; k++) {
				ndxs[k]=dimFragments[mebfFragments[mebf][k][0]][2];
				lens[k]=lenStateList[ndxs[k]];
				if(spinStateList[ndxs[k]]==1 || spinStateList[ndxs[k]]==3) {
					for(int l=0; l<5; l++) sw[k][l]=sw1[l];
				} else {
					for(int l=0; l<5; l++) sw[k][l]=sw2[l];
				}
			}
			if(curr<prev) {
				mebfSpecification[mebf][3]=curr;
			} else {
				count=0;
				for(int j0=0; j0<lens[0]; j0++) {
					for(int j1=0; j1<lens[1]; j1++) {
						for(int j2=0; j2<lens[2]; j2++) {
							include=false;
							int nS=0;
							int nT=0;
							int nD=0;
							int nS1=0;
							int nq=0;
							int nQ=0;
							
							if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("S")) nS++;
							if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("S")) nS++;
							if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("S")) nS++;
							
							if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("T")) nT++;
							if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("T")) nT++;
							if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("T")) nT++;
							
							if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("D")) nD++;
							if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("D")) nD++;
							if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("D")) nD++;
							
							if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("S1")) nS1++;
							if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("S1")) nS1++;
							if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("S1")) nS1++;
							
							if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("q")) nq++;
							if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("q")) nq++;
							if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("q")) nq++;
							
							if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("Q")) nQ++;
							if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("Q")) nQ++;
							if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("Q")) nQ++;
							
							if(nS==3) include=true;
							if(nS==1 && nT==2) include=true;
							if(nS==1 && nD==2) include=true;
							if(nS==1 && nQ==2) include=true;
							
							if(expansion==0 && nD>0) include=false;
							if(expansion<2 && nS1>1) include=false;
							if(expansion<2 && nS1>0 && nD>0) include=false;
							if(expansion<2 && nS1>0 && nT>0) include=false;

							if(charge==0) {
								int chg=0;
								if(stateNames[ndxStateList[ndxs[0]][j0]].trim().indexOf("+")>=0) chg++;
								if(stateNames[ndxStateList[ndxs[1]][j1]].trim().indexOf("+")>=0) chg++;
								if(stateNames[ndxStateList[ndxs[2]][j2]].trim().indexOf("+")>=0) chg++;
								if(stateNames[ndxStateList[ndxs[0]][j0]].trim().indexOf("-")>=0) chg--;
								if(stateNames[ndxStateList[ndxs[1]][j1]].trim().indexOf("-")>=0) chg--;
								if(stateNames[ndxStateList[ndxs[2]][j2]].trim().indexOf("-")>=0) chg--;
								if(chg!=charge) include=false;
							}
							if(include) {
								mebfFragments[mebf][0][count+1]=ndxStateList[ndxs[0]][j0];
								mebfFragments[mebf][1][count+1]=ndxStateList[ndxs[1]][j1];
								mebfFragments[mebf][2][count+1]=ndxStateList[ndxs[2]][j2];
								count++;
							}
						}
					}
				}
				mebfSpecification[mebf][3]=count;
			}
		}

		// For trimer with spin 2 Disabled for now
		if(mebfSpecification[mebf][0]==3 && mebfSpecification[mebf][1]==22) {
			for(int k=0; k<3; k++) {
				ndxs[k]=dimFragments[mebfFragments[mebf][k][0]][2];
				lens[k]=lenStateList[ndxs[k]];
				if(spinStateList[ndxs[k]]==1 || spinStateList[ndxs[k]]==3) {
					for(int l=0; l<3; l++) sw[k][l]=sw1[l];
				} else {
					for(int l=0; l<3; l++) sw[k][l]=sw2[l];
				}
			}
			if(curr<prev) {
				mebfSpecification[mebf][3]=curr;
			} else {
				count=0;
				for(int j0=0; j0<lens[0]; j0++) {
					for(int j1=0; j1<lens[1]; j1++) {
						for(int j2=0; j2<lens[2]; j2++) {
							include=false;
							int nS=0;
							int nT=0;
							int nD=0;
							int nS1=0;
							
							if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("S")) nS++;
							if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("S")) nS++;
							if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("S")) nS++;
							
							if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("T")) nT++;
							if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("T")) nT++;
							if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("T")) nT++;
							
							if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("D")) nD++;
							if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("D")) nD++;
							if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("D")) nD++;
							
							if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("S1")) nS1++;
							if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("S1")) nS1++;
							if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("S1")) nS1++;
							
							if(nS==2 && nD==1) include=true;
							if(nT==2 && nD==1) include=true;
							
							if(expansion==0 && nT>1 &&nD>1) include=false;
							if(expansion<2 && nS1>0 && nT>0) include=false;
							
							if(charge==0) {
								int chg=0;
								if(stateNames[ndxStateList[ndxs[0]][j0]].trim().indexOf("+")>=0) chg++;
								if(stateNames[ndxStateList[ndxs[1]][j1]].trim().indexOf("+")>=0) chg++;
								if(stateNames[ndxStateList[ndxs[2]][j2]].trim().indexOf("+")>=0) chg++;
								if(stateNames[ndxStateList[ndxs[0]][j0]].trim().indexOf("-")>=0) chg--;
								if(stateNames[ndxStateList[ndxs[1]][j1]].trim().indexOf("-")>=0) chg--;
								if(stateNames[ndxStateList[ndxs[2]][j2]].trim().indexOf("-")>=0) chg--;
								if(chg!=charge) include=false;
							}
							if(include) {
								mebfFragments[mebf][0][count+1]=ndxStateList[ndxs[0]][j0];
								mebfFragments[mebf][1][count+1]=ndxStateList[ndxs[1]][j1];
								mebfFragments[mebf][2][count+1]=ndxStateList[ndxs[2]][j2];
								count++;
							}
						}
					}
				}
				mebfSpecification[mebf][3]=count;
			}
		}

		// For trimer with spin 3
		if(mebfSpecification[mebf][0]==3 && mebfSpecification[mebf][1]==3) {
			for(int k=0; k<3; k++) {
				ndxs[k]=dimFragments[mebfFragments[mebf][k][0]][2];
				lens[k]=lenStateList[ndxs[k]];
				if(spinStateList[ndxs[k]]==1 || spinStateList[ndxs[k]]==3) {
					for(int l=0; l<3; l++) sw[k][l]=sw1[l];
				} else {
					for(int l=0; l<3; l++) sw[k][l]=sw2[l];
				}
			}
			if(curr<prev) {
				mebfSpecification[mebf][3]=curr;
			} else {
				count=0;
				for(int j0=0; j0<lens[0]; j0++) {
					for(int j1=0; j1<lens[1]; j1++) {
						for(int j2=0; j2<lens[2]; j2++) {
							include=false;
							int nS=0;
							int nT=0;
							int nD=0;
							int nS1=0;
							
							if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("S")) nS++;
							if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("S")) nS++;
							if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("S")) nS++;
							
							if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("T")) nT++;
							if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("T")) nT++;
							if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("T")) nT++;
							
							if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("D")) nD++;
							if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("D")) nD++;
							if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("D")) nD++;
							
							if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("S1")) nS1++;
							if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("S1")) nS1++;
							if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("S1")) nS1++;
							
							if(nS==2 && nT==1) include=true;
							if(nT==3) include=true;
							if(nT==1 && nD==2) include=true;
							
							if(nS1>1) include=false;
							if(nT>0 && nS1>0) include=false;
							if(nD>0 && nS1>0) include=false;
							
							if(charge==0) {
								int chg=0;
								if(stateNames[ndxStateList[ndxs[0]][j0]].trim().indexOf("+")>=0) chg++;
								if(stateNames[ndxStateList[ndxs[1]][j1]].trim().indexOf("+")>=0) chg++;
								if(stateNames[ndxStateList[ndxs[2]][j2]].trim().indexOf("+")>=0) chg++;
								if(stateNames[ndxStateList[ndxs[0]][j0]].trim().indexOf("-")>=0) chg--;
								if(stateNames[ndxStateList[ndxs[1]][j1]].trim().indexOf("-")>=0) chg--;
								if(stateNames[ndxStateList[ndxs[2]][j2]].trim().indexOf("-")>=0) chg--;
								if(chg!=charge) include=false;
							}
							if(include) {
								mebfFragments[mebf][0][count+1]=ndxStateList[ndxs[0]][j0];
								mebfFragments[mebf][1][count+1]=ndxStateList[ndxs[1]][j1];
								mebfFragments[mebf][2][count+1]=ndxStateList[ndxs[2]][j2];
								count++;
							}
						}
					}
				}
				mebfSpecification[mebf][3]=count;
			}
		}
		
		// For tetramer with spin 1
		if(mebfSpecification[mebf][0]==4 && mebfSpecification[mebf][1]==1) {
			for(int k=0; k<4; k++) {
				ndxs[k]=dimFragments[mebfFragments[mebf][k][0]][2];
				lens[k]=lenStateList[ndxs[k]];
				if(spinStateList[ndxs[k]]==1 || spinStateList[ndxs[k]]==3) {
					for(int l=0; l<3; l++) sw[k][l]=sw1[l];
				} else {
					for(int l=0; l<3; l++) sw[k][l]=sw2[l];
				}
			}
			if(curr<prev) {
				mebfSpecification[mebf][3]=curr;
			} else {
				count=0;
				for(int j0=0; j0<lens[0]; j0++) {
					for(int j1=0; j1<lens[1]; j1++) {
						for(int j2=0; j2<lens[2]; j2++) {
							for(int j3=0; j3<lens[3]; j3++) {
								include=false;
								int nS=0;
								int nT=0;
								int nD=0;
								int nS1=0;
								
								if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("S")) nS++;
								if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("S")) nS++;
								if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("S")) nS++;
								if(stateNames[ndxStateList[ndxs[3]][j3]].trim().startsWith("S")) nS++;
								
								if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("T")) nT++;
								if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("T")) nT++;
								if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("T")) nT++;
								if(stateNames[ndxStateList[ndxs[3]][j3]].trim().startsWith("T")) nT++;
								
								if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("D")) nD++;
								if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("D")) nD++;
								if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("D")) nD++;
								if(stateNames[ndxStateList[ndxs[3]][j3]].trim().startsWith("D")) nD++;
								
								if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("S1")) nS1++;
								if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("S1")) nS1++;
								if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("S1")) nS1++;
								if(stateNames[ndxStateList[ndxs[3]][j3]].trim().startsWith("S1")) nS1++;
								
								if(nS==4) include=true;
								if(nS==2 && nD==2) include=true;
								if(nS==2 && nT==2) include=true;
								if(nD==2 && nT==2) include=true;
								
								if(expansion==0 && nS1>0 && nT>0) include=false;
								if(expansion==0 && nD>0) include=false;
																
								if(expansion<2 && nS1>1) include=false;
								if(expansion<2 && nD>2) include=false;
								if(expansion<2 && nD>0 && nS1>0) include=false;
								if(expansion<2 && nD>0 && nT>0) include=false;
								if(expansion<2 && nS1>0 && nT>0) include=false;
							
								if(charge==0) {
									int chg=0;
									if(stateNames[ndxStateList[ndxs[0]][j0]].trim().indexOf("+")>=0) chg++;
									if(stateNames[ndxStateList[ndxs[1]][j1]].trim().indexOf("+")>=0) chg++;
									if(stateNames[ndxStateList[ndxs[2]][j2]].trim().indexOf("+")>=0) chg++;
									if(stateNames[ndxStateList[ndxs[3]][j3]].trim().indexOf("+")>=0) chg++;
									if(stateNames[ndxStateList[ndxs[0]][j0]].trim().indexOf("-")>=0) chg--;
									if(stateNames[ndxStateList[ndxs[1]][j1]].trim().indexOf("-")>=0) chg--;
									if(stateNames[ndxStateList[ndxs[2]][j2]].trim().indexOf("-")>=0) chg--;
									if(stateNames[ndxStateList[ndxs[3]][j3]].trim().indexOf("-")>=0) chg--;
									if(chg!=charge) include=false;
								}
								if(include && count<maxMEBFstates-1) {
									mebfFragments[mebf][0][count+1]=ndxStateList[ndxs[0]][j0];
									mebfFragments[mebf][1][count+1]=ndxStateList[ndxs[1]][j1];
									mebfFragments[mebf][2][count+1]=ndxStateList[ndxs[2]][j2];
									mebfFragments[mebf][3][count+1]=ndxStateList[ndxs[3]][j3];
									count++;
								}
							}
						}
					}
				}
				mebfSpecification[mebf][3]=count;
			}
		}

		// For pentamer with spin 1
		if(mebfSpecification[mebf][0]==5 && mebfSpecification[mebf][1]==1) {
			for(int k=0; k<5; k++) {
				ndxs[k]=dimFragments[mebfFragments[mebf][k][0]][2];
				lens[k]=lenStateList[ndxs[k]];
				if(spinStateList[ndxs[k]]==1 || spinStateList[ndxs[k]]==3) {
					for(int l=0; l<3; l++) sw[k][l]=sw1[l];
				} else {
					for(int l=0; l<3; l++) sw[k][l]=sw2[l];
				}
			}
			if(curr<prev) {
				mebfSpecification[mebf][3]=curr;
			} else {
				count=0;
				for(int j0=0; j0<lens[0]; j0++) {
					for(int j1=0; j1<lens[1]; j1++) {
						for(int j2=0; j2<lens[2]; j2++) {
							for(int j3=0; j3<lens[3]; j3++) {
								for(int j4=0; j4<lens[4]; j4++) {
									include=false;
									int nS=0;
									int nT=0;
									int nD=0;
									int nS1=0;
									
									if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("S")) nS++;
									if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("S")) nS++;
									if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("S")) nS++;
									if(stateNames[ndxStateList[ndxs[3]][j3]].trim().startsWith("S")) nS++;
									if(stateNames[ndxStateList[ndxs[4]][j4]].trim().startsWith("S")) nS++;
									
									if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("T")) nT++;
									if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("T")) nT++;
									if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("T")) nT++;
									if(stateNames[ndxStateList[ndxs[3]][j3]].trim().startsWith("T")) nT++;
									if(stateNames[ndxStateList[ndxs[4]][j4]].trim().startsWith("T")) nT++;
									
									if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("D")) nD++;
									if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("D")) nD++;
									if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("D")) nD++;
									if(stateNames[ndxStateList[ndxs[3]][j3]].trim().startsWith("D")) nD++;
									if(stateNames[ndxStateList[ndxs[4]][j4]].trim().startsWith("D")) nD++;
									
									if(stateNames[ndxStateList[ndxs[0]][j0]].trim().startsWith("S1")) nS1++;
									if(stateNames[ndxStateList[ndxs[1]][j1]].trim().startsWith("S1")) nS1++;
									if(stateNames[ndxStateList[ndxs[2]][j2]].trim().startsWith("S1")) nS1++;
									if(stateNames[ndxStateList[ndxs[3]][j3]].trim().startsWith("S1")) nS1++;
									if(stateNames[ndxStateList[ndxs[4]][j4]].trim().startsWith("S1")) nS1++;
									
									if(nS==5 && nS1==0) include=true;
									if(nS==5 && nS1==1) include=true;
									if(nS==3 && nD==2) include=true;
									if(nS==2 && nT==2 && nS1==0) include=true;
																	
									if(charge==0) {
										int chg=0;
										if(stateNames[ndxStateList[ndxs[0]][j0]].trim().indexOf("+")>=0) chg++;
										if(stateNames[ndxStateList[ndxs[1]][j1]].trim().indexOf("+")>=0) chg++;
										if(stateNames[ndxStateList[ndxs[2]][j2]].trim().indexOf("+")>=0) chg++;
										if(stateNames[ndxStateList[ndxs[3]][j3]].trim().indexOf("+")>=0) chg++;
										if(stateNames[ndxStateList[ndxs[4]][j4]].trim().indexOf("+")>=0) chg++;
										if(stateNames[ndxStateList[ndxs[0]][j0]].trim().indexOf("-")>=0) chg--;
										if(stateNames[ndxStateList[ndxs[1]][j1]].trim().indexOf("-")>=0) chg--;
										if(stateNames[ndxStateList[ndxs[2]][j2]].trim().indexOf("-")>=0) chg--;
										if(stateNames[ndxStateList[ndxs[3]][j3]].trim().indexOf("-")>=0) chg--;
										if(stateNames[ndxStateList[ndxs[4]][j4]].trim().indexOf("-")>=0) chg--;
										if(chg!=charge) include=false;
									}
									if(include && count<maxMEBFstates-1) {
										mebfFragments[mebf][0][count+1]=ndxStateList[ndxs[0]][j0];
										mebfFragments[mebf][1][count+1]=ndxStateList[ndxs[1]][j1];
										mebfFragments[mebf][2][count+1]=ndxStateList[ndxs[2]][j2];
										mebfFragments[mebf][3][count+1]=ndxStateList[ndxs[3]][j3];
										mebfFragments[mebf][4][count+1]=ndxStateList[ndxs[4]][j4];
										count++;
									}
								}
							}
						}
					}
				}
				mebfSpecification[mebf][3]=count;
			}
		}

		update();
	}
	
	private void mebfCellSelected(MouseEvent e) {
		Integer row = mebfsTable.getSelectedRow();
		Integer col = mebfsTable.getSelectedColumn();
		JFrame jf = new JFrame();
		String value;
		Integer mebf = mebfIndex[row][0];
		Integer nmer = mebfSpecification[mebf][0];
		Integer spin = mebfSpecification[mebf][1];
		Integer chrg = mebfSpecification[mebf][2];
		Integer stat = mebfSpecification[mebf][3];
		if(col==1) {
			try {
				value = JOptionPane.showInputDialog(jf,"Enter number of fragments for MEBF "+mebfName[mebf].trim());
				if(value.length()>0) {
					mebfSpecification[mebf][0]=Integer.valueOf(value);
					if(mebfSpecification[mebf][0]<1) mebfSpecification[mebf][0]=1;
				}
			} catch (NullPointerException e1) {
			}

			selectMEBFStates(mebf,nmer,spin,chrg,stat);
			update();			
		}
		
		if(col==2) {
			try {
				value = JOptionPane.showInputDialog(jf,"Enter spin for MEBF "+mebfName[mebf].trim());
				if(value.length()>0) mebfSpecification[mebf][1]=Integer.valueOf(value);
			} catch (NullPointerException e1) {
			}

			selectMEBFStates(mebf,nmer,spin,chrg, stat);
			update();
		}

		if(col==3) {
			try {
				value = JOptionPane.showInputDialog(jf,"Enter charge for MEBF "+mebfName[mebf].trim());
				if(value.length()>0) mebfSpecification[mebf][2]=Integer.valueOf(value);
			} catch (NullPointerException e1) {
			}

			selectMEBFStates(mebf,nmer,spin,chrg, stat);
			update();
		}
		
		if(col==4) {
			Integer prev=mebfSpecification[mebf][3];
			try {
				value = JOptionPane.showInputDialog(jf,"Enter number of states for MEBF "+mebfName[mebf].trim());
				if(value.length()>0) mebfSpecification[mebf][3]=Integer.valueOf(value);
			} catch (NullPointerException e1) {
		}

			selectMEBFStates(mebf,nmer,spin,chrg, stat);
			update();
		}
	}
	
	private void mebfEnergiesCellSelected(MouseEvent e) {
	}

	private void nociEnergiesCellSelected(MouseEvent e) {
	}
	
	private void write_Molcas_MEBF_files() {
		for(int i=0; i<numMEBFs; i++) {
			Integer nfrags = mebfSpecification[i][0];
			String pName = " ";
			String[] frags = new String[nfrags];
			String[] source = new String[nfrags];
			String[] fname = new String[nfrags];
			Integer[] fstat = new Integer[nfrags];
			Integer[] nums = new Integer[nfrags];
			Double[][] randt = new Double[nfrags][6];
			String fileName = projectName.trim()+mebfName[i].trim();
			if(nfrags==1) fileName=fileName.trim()+mebfName[i].trim();
			for(int j=0; j<nfrags; j++) {
				int last=namFragments[j].indexOf("_");
				fname[j]=namFragments[j].substring(0,last);
				frags[j]=fragmentNames[mebfFragments[i][j][0]];
				source[j]=fragmentNames[mebfFragments[i][j][1]];
//				System.out.println("FRAG "+j+" : "+frags[j]+" SOURCE "+source[j]);
				fstat[j]=dimFragments[j][2];
				nums[j]=lenStateList[dimFragments[j][2]];
				pName = projectName.trim();
				for(int k=0; k<6; k++) randt[j][k]=movFragments[mebfFragments[i][j][0]][k];
			}
			if(fragment.write_MEBF_XYZ(fileName, pName, nfrags, fname, frags, randt)) {
				fragment.write_Molcas_MEBF_One(fileName, pName, nfrags, fname, frags, randt, basisSets[basisSet], contracts[contract], cholesky);	
				fragment.write_Molcas_MEBF_CB(fileName,projectName, nfrags, frags, fstat, lenStateList, ndxStateList, thresh_MO);	
				fragment.write_Molcas_MEBF_Two(fileName, pName, nfrags, fname, frags, randt, basisSets[basisSet], contracts[contract], cholesky);
				fragment.write_Run_Script_MEBFs(fileName, pName, nfrags, fname,frags,source,fstat, lenStateList, ndxStateList, numRanks,memory,fragmentDefinitions,account,jobName,timeLimit);
			}
		}
	}

	public void write_GronOR_NOCI() {
		Integer numME =0;
		Integer nmer = 0;
		Integer spin = 1;
		Integer isp1, isp2, isp3, isp4;
		
		for(int i=0; i<5; i++) {
			for(int j=0; j<5; j++) {
				for(int k=0; k<5; k++) {
					for(int m=0; m<5; m++) {
						couple3[i][j][k][m]=9;
					}
				}
			}
		}
		couple3[0][0][0][0]=1; // S,S,S
		couple3[0][0][1][1]=2; // S,D,D
		couple3[0][1][0][1]=2; // D,S,D
		couple3[0][1][1][0]=1; // D,D,S
		couple3[0][0][2][2]=3; // S,T,T
		couple3[0][2][0][2]=3; // T,S,T
		couple3[0][2][2][0]=1; // T,T,S
		couple3[0][0][4][4]=5; // S,Q,Q
		couple3[0][4][0][4]=5; // Q,S,Q
		couple3[0][4][4][0]=1; // Q,Q,S
		
		couple3[1][1][0][0]=2; // D,S,S
		couple3[1][0][1][0]=2; // S,D,S
		couple3[1][0][0][1]=1; // S,S,D
		couple3[1][1][1][1]=1; // D,D,D
		
		couple3[2][0][0][2]=1; // S,S,T
		couple3[2][0][2][0]=3; // S,T,S
		couple3[2][2][0][0]=3; // T,S,S
		couple3[2][2][2][2]=1; // T,T,T
		
		couple3[2][1][1][2]=1; // D,D,T
		couple3[2][1][2][1]=2; // D,T,D
		couple3[2][2][1][1]=2; // T,D,D
		
		couple3[2][2][2][2]=1; // T,T,T

		for(int i=0; i<5; i++) {
			for(int j=0; j<5; j++) {
				for(int k=0; k<5; k++) {
					for(int m=0; m<5; m++) {
						for(int n=0; n<5; n++) {
							for(int l=0; l<2; l++) {
								couple4[i][j][k][m][n][l]=9;
							}
						}
					}
				}
			}
		}
		couple4[0][0][0][0][0][0]=1; // S,S,S,S
		couple4[0][0][0][0][0][1]=1;
		
		couple4[0][0][0][1][1][0]=1; // S,S,D,D
		couple4[0][0][0][1][1][1]=2;
		couple4[0][0][1][0][1][0]=2; // S,D,S,D
		couple4[0][0][1][0][1][1]=2;
		couple4[0][0][1][1][0][0]=2; // S,D,D,S
		couple4[0][0][1][1][0][1]=1;
		couple4[0][1][0][0][1][0]=2; // D,S,S,D
		couple4[0][1][0][0][1][1]=2;
		couple4[0][1][0][1][0][0]=2; // D,S,D,S
		couple4[0][1][0][1][0][1]=1;
		couple4[0][1][1][0][0][0]=1; // D,D,S,S
		couple4[0][1][1][0][0][1]=1;
		
		couple4[0][0][0][2][2][0]=1; // S,S,T,T
		couple4[0][0][0][2][2][1]=3;
		couple4[0][0][2][0][2][0]=3; // S,T,S,T
		couple4[0][0][2][0][2][1]=3;
		couple4[0][0][2][2][0][0]=3; // S,T,T,S
		couple4[0][0][2][2][0][1]=1;
		couple4[0][2][0][0][2][0]=3; // T,S,S,T
		couple4[0][2][0][0][2][1]=3;
		couple4[0][2][0][2][0][0]=3; // T,S,T,S
		couple4[0][2][0][2][0][1]=1;
		couple4[0][2][2][0][0][0]=1; // T,T,S,S
		couple4[0][2][2][0][0][1]=1;
		
		
		for(int i=0; i<numMEBFs; i++) {
			String fileName = projectName.trim()+mebfName[i].trim()+"_GronOR.inp";
			numME=mebfSpecification[i][3];
			nmer=mebfSpecification[i][0];
			spin=mebfSpecification[i][1];
			try {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("MEBFs "+projectName.trim()+" "+numME);
				for(int j=0; j<nmer; j++) {
					inputFile.print(fragmentNames[mebfFragments[i][j][0]]+" ");
					for(int k=0; k<numME; k++) inputFile.print(" "+stateNames[mebfFragments[i][j][k+1]]);
					inputFile.println();
				}
				if(nmer==3) {
					inputFile.println("Couplings ");
					inputFile.print("  ");
					for(int k=0; k<numME; k++) {
						isp1=stateSpins[mebfFragments[i][0][k+1]]-1;
						isp2=stateSpins[mebfFragments[i][1][k+1]]-1;
						isp3=stateSpins[mebfFragments[i][2][k+1]]-1;
						inputFile.print("  "+couple3[spin-1][isp1][isp2][isp3]);
					}
					inputFile.println();
				}
				if(nmer==4) {
					inputFile.println("Couplings ");
					for(int l=0; l<2; l++) {
						inputFile.print("  ");
						for(int k=0; k<numME; k++) {
							isp1=stateSpins[mebfFragments[i][0][k+1]]-1;
							isp2=stateSpins[mebfFragments[i][1][k+1]]-1;
							isp3=stateSpins[mebfFragments[i][2][k+1]]-1;
							isp4=stateSpins[mebfFragments[i][3][k+1]]-1;
							inputFile.print("  "+couple4[spin-1][isp1][isp2][isp3][isp4][l]);
						}
						inputFile.println();
					}
				}
//				inputFile.println("Ecorr");
//				for(int j=0; j<nmer; j++) {
////					m=dimFragments[j][0];
//					int mf=mebfFragments[i][j][0];
//					int ms=dimFragments[mf][2];
//					int ls=lenStateList[ms];
//					for(int k=0; k<ls; k++) inputFile.print(" "+(energiesCASPT2[mf][k]-energiesCASSCF[mf][k]));
//					inputFile.println();
//				}
				inputFile.println("Spin "+spin);
				inputFile.println("Threshold "+thresh_CI);
				inputFile.println("Print medium");
				inputFile.close();
			} catch(IOException e) {
			}
		}
	}
	
	public Boolean read_GronOR_arx(Integer mebf) {
		Integer numE=0;
		Integer numC=0;
		Double eMEBF=0.0;
		Double eNOCI=0.0;
		Integer ndx=0;
		numEnergiesGronOR[mebf][0]=0;
		numEnergiesGronOR[mebf][1]=0;
		String fileName = projectName.trim()+mebfName[mebf].trim()+"_GronOR.arx";
		try {
	    	BufferedReader br = new BufferedReader(new FileReader(fileName));
	    	String card;
	    	while((card=br.readLine()) != null) {
				if(card.startsWith("Hamil")) {
					numE=Integer.valueOf(card.substring(5,15).trim());
					numC=Integer.valueOf(card.substring(15,25).trim());
					card=br.readLine();
					if(numE<=numC) {
						for(int j=0; j<numE; j++) {
							card=br.readLine();
							eMEBF=Double.valueOf(card.substring(20*j+6,20*j+26)).doubleValue();
							energiesGronOR[j][mebf][0]=eMEBF;
							numEnergiesGronOR[mebf][0]=j+1;
						}
					} else if(numE<=2*numC) {
						for(int j=0; j<numE; j++) {
							card=br.readLine();
							if(j>=numC) card=br.readLine();
							ndx=j;
							if(ndx>=numC) ndx=numC;
							eMEBF=Double.valueOf(card.substring(20*ndx+6,20*ndx+26)).doubleValue();
							energiesGronOR[j][mebf][0]=eMEBF;
							numEnergiesGronOR[mebf][0]=j+1;
							if(j<numC) card=br.readLine();
						}
					}
				}
				if(card.startsWith("NOCI")) {
					numE=Integer.valueOf(card.substring(5,15).trim());
					numC=Integer.valueOf(card.substring(15,25).trim());
					card=br.readLine();
					if(numE<=numC) {
						card=br.readLine();
						for(int j=0; j<numE; j++) {
							eNOCI=Double.valueOf(card.substring(20*j,20*j+20)).doubleValue();
							energiesGronOR[j][mebf][1]=eNOCI;
							numEnergiesGronOR[mebf][1]=j+1;
						}
						card=br.readLine();
						for(int j=0; j<numE; j++) {
							card=br.readLine();
							for(int k=0; k<numE; k++) {
								coefGronOR[j][k][mebf]=Double.valueOf(card.substring(20*k,20*k+20)).doubleValue();
							}
						}
					} else if(numE<=2*numC) {
						card=br.readLine();
						for(int j=0; j<7; j++) {
							eNOCI=Double.valueOf(card.substring(20*j,20*j+20)).doubleValue();
							energiesGronOR[j][mebf][1]=eNOCI;
							numEnergiesGronOR[mebf][1]=j+1;
						}
						card=br.readLine();
						for(int j=0; j<numE; j++) {
							card=br.readLine();
							for(int k=0; k<7; k++) {
								coefGronOR[j][k][mebf]=Double.valueOf(card.substring(20*k,20*k+20)).doubleValue();
							}
						}
						card=br.readLine();
						for(int j=7; j<numE; j++) {
							ndx=j-7;
							eNOCI=Double.valueOf(card.substring(20*ndx,20*ndx+20)).doubleValue();
							energiesGronOR[j][mebf][1]=eNOCI;
							numEnergiesGronOR[mebf][1]=j+1;
						}
						card=br.readLine();
						for(int j=0; j<numE; j++) {
							card=br.readLine();
							for(int k=7; k<numE; k++) {
								ndx=k-7;
								coefGronOR[j][k][mebf]=Double.valueOf(card.substring(20*ndx,20*ndx+20)).doubleValue();
							}
						}
					}
				}
	    	}
	    	br.close();
		} catch(Exception ee) {
			return false;
		}
    	return true;
	}
	
	public void actionPerformed(ActionEvent e){}

	public void windowClosing(WindowEvent event){ setVisible(false); System.exit(0);}

	public void windowClosed(WindowEvent event) { System.exit(0); }
	  
	public void windowDeiconified(WindowEvent event) {}
	  
	public void windowIconified(WindowEvent event) {}
	  
	public void windowActivated(WindowEvent event) {}
	  
	public void windowDeactivated(WindowEvent event) {}
	  
	public void windowOpened(WindowEvent event) {}	

    public void stateChanged(ChangeEvent e) {}

    public void mouseClicked(MouseEvent mouse) {}

    public void mousePressed(MouseEvent mouse){}

    public void mouseReleased(MouseEvent mouse){}

    public void mouseEntered(MouseEvent mouse){}

    public void mouseExited(MouseEvent mouse){}
    
    public void itemStateChanged(ItemEvent item){}
    
}
