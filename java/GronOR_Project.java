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

public class GronOR_Project extends JFrame implements ActionListener, ChangeListener, WindowListener, MouseListener{
	
	private static final long serialVersionUID = 4L;

	JMenuBar menubar = new JMenuBar();
	JMenuItem create, open, close, list, save, quit;
	String projectName, projectFile;

	JPanel parametersPanel;
	JPanel dimensionPanel;
	JPanel numberPanel;
	JPanel threshPanel;
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
	Integer[] fragmentStates = new Integer[25];					// Index to stateNames[] for each unique fragment state used
	
	Integer[] ndxMebfTable = new Integer[25];					// Index of each state into into MEBF table 
	Integer[] ndxMebfState = new Integer[25];					// Index into MEBF table for each state in MEBF 
	
	Integer[][] dimFragments = new Integer[maxFragments][11];	// Dimensions of fragments: number of atoms, states, electrons, etc.
	String[] namFragments = new String[maxFragments];			// Names of fragments used to determine xyz-formatted coordinate origin 
	Double[][] movFragments = new Double[maxFragments][6];		// Rotation and translation of coordinates with respect to original source	
	Double[] energiesDFT = new Double[maxFragments];			// DFT optimized energy of S0 state from NWChem
	Double[] energiesSCF = new Double[maxFragments];			// SCF energy of S0 state from Molcas
	Double[][] energiesCASSCF = new Double[maxFragments][12];	// CASSCF energies of all states of fragment from Molcas
	Double[][] energiesCASPT2 = new Double[maxFragments][12];	// CASPT2 energies of all states of fragment from Molcas
	
	
	Integer numRanks=12;		// Number of ranks in internal mpirun
	Integer numTBD1=0;
	Integer numTBD2=0;
	Integer numTBD3=0;
	
	Integer numEnergies = 0;	// Number of energy entries

	Integer numMEBFRows = 0;
    Integer numMEBFs = 0;
    Integer newMEBFs = 0;
    Integer maxMEBFs = 32;
    Integer maxMer = 3;
    Integer maxMEBFstates=18;
    
    Integer numStateLabels = 12;

	Integer[][] numEnergiesGronOR = new Integer[maxMEBFs][2];		// Number of GronORenergy entries
	Double[][][] energiesGronOR = new Double[25][maxMEBFs][2];
	Double[][][] energiesRelGronOR = new Double[25][maxMEBFs][2];
	Double[][][] coefGronOR = new Double[25][25][maxMEBFs];

    String[] stateListLabels = new String[] {"ID","Num","Spin","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"};
    String[] fragmentLabels = new String[] {"ID", "SRC", "XYZ File", "Atoms", "States", "Electrons", "CASe", "CASo", "Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "Alt"};
    String[] fragmentNames = new String[] {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"};
    String[] energyNames = new String[] {"E(DFT)", "E(SCF)", "E(CASSCF)", "E(CASPT2)"};
    String[] stateLabels = new String[] {" ", "S0", "S1", "T1", "D-", "D+", "S2", "T2", "E-", "E+"," "," "};

    String[] stateNames = new String[] {"S0","S1","S2","D0","D1","T1","T2","S+","D+","T+","S-","D-","T-"};
	String[] mebfLabels = new String[] {"ID","n-mer","Spin","Charge","States","Frag","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19" };
	String[] nociLabels = new String[] {"ID", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9" };

	Integer[] stateLabelIndex = new Integer[numStateLabels+3];
	
    JTable statesTable;
    JTable fragmentsTable;
    JTable energiesTable;
    JTable mebfsTable;
	JTable dimensionTable = new JTable();
	JTable numberTable = new JTable();
	JTable threshTable = new JTable();
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
    Object[][] mebfDefinitions = new Object[maxMEBFs*maxMer][25];
    Object[][] mebfEnergies = new Object[maxMEBFs][10];
    Object[][] nociEnergies = new Object[25][maxMEBFs];
    
    Integer[][] mebfIndex = new Integer[maxMEBFs*maxMer][2];
    Integer[][] mebfSpecification = new Integer[maxMEBFs][5];	 		// 0: n-mer; 1: spin; 2: charge; 3:number of states to include 4:fragment; 
    String[] mebfName = new String[maxMEBFs];							// name of mefb, e.g. AB for fragment A,B combination
    Integer[][][] mebfFragments = new Integer[maxMEBFs][maxMer][19];	// index to fragment states included in this mebf
    
    Integer numberStateEnergies = 0;
    Integer numberMebfDefinitions = 0;
    Integer numberMebfEnergies = 0;
    Integer numberFragmentDefinitions = 0;
    Integer numberStateDefinitions = 0;
    
    Integer[][] energyFragment = new Integer[maxFragments][2];

    GronOR_Fragment fragment = new GronOR_Fragment();

    Double[] RandT = new Double[6];
    
    Double thresh_CI = 1.0e-5;
    Double thresh_MO = 1.0e-5;

    Integer noci0=0;
    Integer noci1=0;
	
	
    JFrame dialogFrame;

    JFileChooser chooser;
    ExtensionFilter filter;
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

	GronOR_Project(String pName){
		
		super("GronOR Project "+pName);
	    super.setSize(1000,900);
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
			clearFile.println("rm -f "+projectName.trim()+"*.input");
			clearFile.println("rm -f "+projectName.trim()+"*.inp");
			clearFile.println("rm -f "+projectName.trim()+"*.nw");
			clearFile.println("rm -f "+projectName.trim()+"*.run");
			clearFile.println("rm -f "+projectName.trim()+"*.xyz");
			clearFile.close();
		} catch(IOException ei) {
		}
		fileName = "clearall.run";
		try {
			PrintfWriter clearFile = new PrintfWriter(new FileWriter(fileName));
			clearFile.println("#!/usr/bin/tcsh");
			clearFile.println("rm -f *.input");
			clearFile.println("rm -f *.inp");
			clearFile.println("rm -f *.nw");
			clearFile.println("rm -f *.LprOrb");
			clearFile.println("rm -f *.ONEINT");
			clearFile.println("rm -f *.RasOrb*");
			clearFile.println("rm -f *.INPORB*");
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
			clearFile.println("rm -f *.sym");
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
			clearFile.println("rm -f *A.xyz");
			clearFile.println("rm -f *B.xyz");
			clearFile.println("rm -f ONEINT*");
			clearFile.println("rm -f RUNFIL*");
			clearFile.println("rm -f TRA*");
			clearFile.println("rm -f xml*");
			clearFile.println("rm -f .project");
			clearFile.println("rm -f *.arx");
			clearFile.println("rm -f *.inp");
			clearFile.println("rm -f *.out");
			clearFile.println("rm -f *.day");
			clearFile.println("rm -f *.cml");
			clearFile.println("rm -f *.pro");
			clearFile.println("rm -f *.rnk");
			clearFile.println("rm -f *.tst");
			clearFile.println("rm -f *.cpr");
			clearFile.println("rm -f *.xmldump");
			clearFile.println("mv "+projectName.trim()+".prj KEEPTHIS");
			clearFile.println("rm -f "+projectName.trim()+"*");
			clearFile.println("mv KEEPTHIS "+projectName.trim()+".prj");
			clearFile.close();
		} catch(IOException ei) {
		}
	}
	
	private Boolean readProjectFile(String fileName) {
		try {
		    	BufferedReader br = new BufferedReader(new FileReader(fileName));
		    	String card;
		    	card=br.readLine();
			    numSets = Integer.valueOf(card.substring(0,6).trim());
		    	newSets = numSets;
			    numFragments = Integer.valueOf(card.substring(6,12).trim());
		    	newFragments = numFragments;
		    	numMEBFs = Integer.valueOf(card.substring(12,18).trim());
		    	newMEBFs = numMEBFs;
		    	numRanks = Integer.valueOf(card.substring(18,24).trim());
		    	thresh_MO=Double.valueOf(card.substring(24,44)).doubleValue();
		    	thresh_CI=Double.valueOf(card.substring(44,64)).doubleValue();
			    
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
						for(int k=0; k<19; k++) mebfFragments[i][j][k]=Integer.valueOf(card.substring(k*6,k*6+6).trim());
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
			    fw.printf("%6d",newSets);
			    fw.printf("%6d",newFragments);
			    fw.printf("%6d",numMEBFs);
			    fw.printf("%6d",numRanks);
			    fw.printf("%20.10f",thresh_MO);
			    fw.printf("%20.10f",thresh_CI);
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
						for(int k=0; k<19; k++) fw.printf("%6d",mebfFragments[i][j][k]); fw.println();
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
	    ExtensionFilter filter;
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
	    	filter = new ExtensionFilter(extension);
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
				if(aL.equals("H")) number=number+1;
				if(aL.equals("C")) number=number+6;
				if(aL.equals("N")) number=number+7;
				if(aL.equals("O")) number=number+8;
				if(aL.equals("F")) number=number+9;
			}
			br.close();
			return number;
		} catch(IOException ef) {
			System.out.println("IOException in XYZ file "+fileName);
			return 0;
		}
	}

	private void update() {

		updateStatesList();
		updateFragmentList();
		updateFragmentDefinitions();
		readOutputFiles();
		updateFragmentEnergies();
		updateMEBFDefinitions();
		updateNOCIEnergies();
		
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
	    String[] stateNames = new String[] {"S0","S1","S2","D0","D1","T1","T2","S+","D+","T+","S-","D-","T-"};
		
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
			fragmentDefinitions[i][1]=fragmentDefinitions[dimFragments[i][0]][0];
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
		
		// Do not repeatedly apply rotations and translation to the same fragment
		for(int i=0; i<numFragments; i++) {
			if(fragmentNames[i].trim().equals(fragmentDefinitions[dimFragments[i][0]][0])) {
				for(int j=0; j<6; j++) movFragments[i][j]=0.0;
			}
		}
		
		for(int i=0; i<numFragments; i++) {
			fragmentDefinitions[i][1]=fragmentDefinitions[dimFragments[i][0]][0];
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
				nameF=nameP.trim();
			}
			fragment.initialize(nameP,nameA,nameF,nameB,dimFragments[i][3],RandT);

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
			
			if(read_GronOR_arx()) {
			
				for(int i=0; i<numMEBFs; i++) {
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
				}
				for(int i=0; i<maxMEBFStates; i++) nociEnergiesTableModel.setValueAt("E"+(i+1),i,0);
				
				
				Double Emin=0.0;
				Double Emax=0.0;
				for(int i=0; i<numMEBFs; i++) {
					for(int j=0; j<numEnergiesGronOR[i][0]; j++) {
						Emin=Math.min(Emin, energiesGronOR[j][i][0]);
					}
					for(int j=0; j<numEnergiesGronOR[i][1]; j++) {
						Emin=Math.min(Emin, energiesGronOR[j][i][1]);
					}
				}
				Emax=Emin;
				for(int i=0; i<numMEBFs; i++) {
					for(int j=0; j<numEnergiesGronOR[i][0]; j++) {
						Emax=Math.max(Emax, energiesGronOR[j][i][0]);
					}
					for(int j=0; j<numEnergiesGronOR[i][1]; j++) {
						Emax=Math.max(Emax, energiesGronOR[j][i][1]);
					}
				}
				for(int i=0; i<numMEBFs; i++) {
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
				for(int i=0; i<numMEBFs; i++) {
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
				for(int i=0; i<numMEBFs; i++) {
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
	
	private void writeInputFiles() {

		String nameA, nameP, nameF, nameS;
		Integer numCASe=0, numCASo=0;
		Boolean withCASPT2=false;
		
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
			if(stateNames[ndxStateList[i][0]].startsWith("D")) mult=2;
			if(stateNames[ndxStateList[i][0]].startsWith("T")) mult=3;
			nameF=projectName.trim()+nameA.trim();
			fragment.write_NWChem_DFT(nameF,mult,numRanks);
		}

		for(int i=0; i<=numRunMolcas; i++) {
			for(int j=0; j<6; j++) RandT[j]=movFragments[i][j];
			nameP=projectName.trim();
			nameA= (String) fragmentDefinitions[i][0];
			stateIndex=dimFragments[i][2];
			withCASPT2=(fragmentDefinitions[i][0]==fragmentDefinitions[i][1]);

			Integer mult=1;
			if(stateNames[ndxStateList[stateIndex][0]].startsWith("D")) mult=2;
			if(stateNames[ndxStateList[stateIndex][0]].startsWith("T")) mult=3;
			nameF=projectName.trim()+nameA.trim();
			fragment.write_Molcas_SCF(nameF,mult);
			
			if(lenStateList[stateIndex]>0) {
				for(int j=0; j<lenStateList[stateIndex]; j++) {
					nameS=stateNames[ndxStateList[stateIndex][j]];
					numCASe=dimFragments[i][4];
					numCASo=dimFragments[i][5];
					fragment.write_Molcas_CASSCF(nameF,nameS,withCASPT2,numCASe,numCASo);

				}		
			}
		}
		
		for(int i=0; i<numFragments; i++) {
			fragment.write_Run_Script_Fragments(i,dimFragments[i][2],lenStateList,ndxStateList,numRanks);
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

			fragment.initialize(nameP, nameA, nameF, nameB, dimFragments[i][3], RandT);
			
			dimFragments[i][6]=0;
			dimFragments[i][7]=0;
			dimFragments[i][8]=0;
			dimFragments[i][9]=0;
			
			if(fragment.NWChem_Converged(i)) {
				energy=fragment.NWChem_DFT(i);
				energiesDFT[i]=energy;
				dimFragments[i][6]=1;
			} 
			
			if(fragment.Molcas_SCF_Converged(i,dimFragments[i][4])) {
				energy=fragment.Molcas_SCF(i,dimFragments[i][4]);
				numAlt=fragment.Molcas_numAlt();
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
		read_GronOR_arx();
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
					for(int k=1; k<19; k++) {
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
			for(int k=0; k<25; k++) {
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
			for(int j=0; j<25; j++) {
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

		dimensionPanel = new JPanel();
		dimensionPanel.setLayout(new BoxLayout(dimensionPanel,BoxLayout.X_AXIS));
		dimensionPanel.setPreferredSize(new Dimension(120,70));
		dimensionPanel.setMinimumSize(new Dimension(120,70));
		dimensionPanel.setMaximumSize(new Dimension(120,70));
		LineBorder dimensionBorder = new LineBorder(Color.black);
		dimensionPanel.setBorder(dimensionBorder);
		Object[][] dimensionData = new Object[][] {
			{"State Sets",numSets},
			{"Fragments",numFragments},
			{"MEBFs",numMEBFs}
		};
		String[] dimensionColumns = new String[] {" "," "};
		dimensionTable = new JTable(dimensionData,dimensionColumns);
		dimensionTable.setCellSelectionEnabled(true);
		dimensionTable.getColumnModel().getColumn(0).setMaxWidth(80);
		dimensionTable.getColumnModel().getColumn(1).setMaxWidth(30);
		ListSelectionModel dimensionSelectionModel = dimensionTable.getSelectionModel();
		dimensionSelectionModel.addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent e) {
				JFrame jf = new JFrame();
				String value;
				if(dimensionTable.getSelectedRow()==0) {
					try {
						value = JOptionPane.showInputDialog(jf,"Enter new number of state definitions");
						if(value.length()>0) newSets=Integer.valueOf(value);
					} catch(NullPointerException e1) {
						newSets=numSets;
					}
					update();
				}
				if(dimensionTable.getSelectedRow()==1) {
					try {
						value = JOptionPane.showInputDialog(jf,"Enter new number of fragments");
						if(value.length()>0) newFragments=Integer.valueOf(value);
					} catch(NullPointerException e1) {
						newFragments=numFragments;
					}
					update();
				}
				if(dimensionTable.getSelectedRow()==2) {
					try {
						value = JOptionPane.showInputDialog(jf,"Enter new number of MEBFs");
						if(value.length()>0) newMEBFs=Integer.valueOf(value);
					} catch(NullPointerException e1) {
						newMEBFs=numMEBFs;
					}
					update();
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
			{"TBD", numTBD1},
			{"TBD", numTBD2}
		};
		String[] numberColumns = new String[] {" "," "};
		numberTable = new JTable(numberData,numberColumns);
		numberTable.setCellSelectionEnabled(true);
		numberTable.getColumnModel().getColumn(0).setMaxWidth(80);
		numberTable.getColumnModel().getColumn(1).setMaxWidth(30);
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
				numberData[0][1]=numRanks;
				numberData[1][1]=0;
				numberData[2][1]=0;
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
			{"TBD", numTBD3}
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
					if(value.length()>0) thresh_MO=Double.valueOf(value);
					} catch(NullPointerException e1) {
					}
				}
				threshData[0][1]=thresh_MO;
				threshData[1][1]=thresh_CI;
				threshData[2][1]=0;
				update();
			}
		});
		
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
			}
		});
		dimensionPanel.add(Box.createRigidArea(new Dimension(5,0)));
		dimensionPanel.add(dimensionTable);
		numberPanel.add(Box.createRigidArea(new Dimension(5,0)));
		numberPanel.add(numberTable);
		threshPanel.add(Box.createRigidArea(new Dimension(5,0)));
		threshPanel.add(threshTable);
		buttonPanel.add(Box.createRigidArea(new Dimension(5,5)));
		buttonPanel.add(updateButton);
		buttonPanel.add(generateButton);
		buttonPanel.add(Box.createVerticalGlue());
		parametersPanel.add(dimensionPanel);
		parametersPanel.add(Box.createRigidArea(new Dimension(5,0)));
		parametersPanel.add(numberPanel);
		parametersPanel.add(Box.createRigidArea(new Dimension(5,0)));
		parametersPanel.add(threshPanel);
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
		setCombo.addMouseListener(new MouseAdapter() {
			public void mouseReleased(MouseEvent e) {
				Integer row = statesTable.getSelectedRow();
				Integer col = statesTable.getSelectedColumn();
				if(row>=0 && col>=3) {
					ndxStateList[row][(col-2)]=setCombo.getSelectedIndex();
				}
				update();
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
		ListSelectionModel fragmentSelectionModel = fragmentsTable.getSelectionModel();
		fragmentSelectionModel.addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent e) {
				Integer row = fragmentsTable.getSelectedRow();
				Integer col = fragmentsTable.getSelectedColumn();
				JFrame jf = new JFrame();
				String value;
				String nameA, nameB, nameP, nameF;
				// change fragment equivalence
				if(col==1) {
					String newID = fragmentDefinitions[row][col].toString().trim();
					for(int i=0; i<numFragments; i++) {
						if(newID.equals(fragmentDefinitions[i][0].toString().trim())){
							dimFragments[row][0]=i;
							namFragments[row]=namFragments[i];
							dimFragments[row][1]=getNumAxyz(namFragments[row]);
							dimFragments[row][3]=getNumExyz(namFragments[row]);
							Integer numXA=dimFragments[row][1]-getNumHxyz(fragmentDefinitions[row][col].toString());
							if(numXA<dimFragments[row][4]) {
								numCASe=Math.max(4,numXA);
								numCASo=Math.max(4,numXA);
								dimFragments[row][4]=numCASe;
								dimFragments[row][5]=numCASo;
							}
						}
					}
				}
				// xyz coordinate file
				if(col==2) {
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
					fragment.initialize(nameP,nameA,nameF,nameB,dimFragments[row][3],RandT);
				}
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
		TitledBorder mebfsBorder = new TitledBorder(new LineBorder(Color.black),"MEBF Definitions");
		mebfsBorder.setTitleColor(Color.black);
		mebfsPanel.setBorder(mebfsBorder);
		mebfsTable = new JTable(mebfsTableModel);
		mebfsTable.getColumnModel().getColumn(0).setMaxWidth(60);
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
		fragmentCombo.addMouseListener(new MouseAdapter() {
			public void mouseReleased(MouseEvent e) {
				Integer row = mebfsTable.getSelectedRow();
				Integer col = mebfsTable.getSelectedColumn();
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
		for(int i=0; i<19; i++) {
			TableColumn stateColumn = mebfsTable.getColumnModel().getColumn(6+i);
			stateColumn.setCellEditor(new DefaultCellEditor(stateCombo));
		}
		stateCombo.addMouseListener(new MouseAdapter() {
			public void mouseReleased(MouseEvent e) {
				Integer row = mebfsTable.getSelectedRow();
				Integer col = mebfsTable.getSelectedColumn();
				if(row>=0 && col>=5) {
					Integer mymebf = mebfIndex[row][0];
					Integer mynmer = mebfIndex[row][1];
					mebfFragments[mymebf][mynmer][col-4] = stateCombo.getSelectedIndex();
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
		Integer[] lens = new Integer[2];
		Integer[] ndxs = new Integer[2];
		String[] sw1 = new String[] {"S","D","T"};
		String[] sw2 = new String[] {"D","S","T"};
		String[][] sw = new String[2][3];
		Boolean include = false;
		
		if(mebfSpecification[mebf][0]>maxMer) mebfSpecification[mebf][0]=maxMer;
		if(mebfSpecification[mebf][0]>nmer) {
			for(int i=nmer; i<mebfSpecification[mebf][0]; i++) {
				for(int k=0; k<19; k++) {
					mebfFragments[mebf][i][k]=0;
				}
				for(int k=0; k<mebfSpecification[mebf][4]+1; k++) {
					mebfFragments[mebf][i][k]=mebfFragments[mebf][nmer-1][k];
				}
				mebfFragments[mebf][i][0]=mebfFragments[mebf][i-1][0]+1;
				if(mebfFragments[mebf][i][0]>maxMer) mebfFragments[mebf][i][0]=mebfFragments[mebf][i][0]-maxMer;
			}
		}
		nmer=mebfSpecification[mebf][0];
		
		spin=mebfSpecification[mebf][1];

		charge=mebfSpecification[mebf][2];

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
				if(spinStateList[ndxs[k]]==1 || spinStateList[ndxs[k]]==3) {
					for(int l=0; l<3; l++) sw[k][l]=sw1[l];
				} else {
					for(int l=0; l<3; l++) sw[k][l]=sw2[l];
				}
			}
			count=0;
			for(int i=0; i<3; i++) {
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
			Integer[] fstat = new Integer[nfrags];
			Integer[] nums = new Integer[nfrags];
			Double[][] randt = new Double[nfrags][6];
			String fileName = projectName.trim()+mebfName[i].trim();
			if(nfrags==1) fileName=fileName.trim()+mebfName[i].trim();
			for(int j=0; j<nfrags; j++) {
				frags[j]=fragmentNames[mebfFragments[i][j][0]];
				fstat[j]=dimFragments[j][2];
				nums[j]=lenStateList[dimFragments[j][2]];
				pName = projectName.trim();
				for(int k=0; k<6; k++) randt[j][k]=movFragments[mebfFragments[i][j][0]][k];
			}
			if(fragment.write_MEBF_XYZ(fileName, pName, nfrags, frags)) {
				fragment.write_Molcas_MEBF_One(fileName, pName, nfrags, frags);	
				fragment.write_Molcas_MEBF_CB(fileName, nfrags, frags, fstat, lenStateList, ndxStateList, thresh_MO);	
				fragment.write_Molcas_MEBF_Two(fileName, pName, nfrags, frags);
				fragment.write_Run_Script_MEBFs(fileName, pName, nfrags, frags, fstat, lenStateList, ndxStateList, numRanks);
			}
		}
	}

	public void write_GronOR_NOCI() {
		Integer numME =0;
		Integer nmer = 0;
		for(int i=0; i<numMEBFs; i++) {
			String fileName = projectName.trim()+mebfName[i].trim()+"_GronOR.inp";
			numME=mebfSpecification[i][3];
			nmer=mebfSpecification[i][0];
			try {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("MEBFs "+projectName.trim()+mebfName[i].trim()+" "+numME);
				for(int j=0; j<nmer; j++) {
					inputFile.print(fragmentNames[mebfFragments[i][j][0]]+" ");
					for(int k=0; k<numME; k++) inputFile.print(" "+stateNames[mebfFragments[i][j][k+1]]);
					inputFile.println();
				}
				inputFile.println("Threshold "+thresh_CI);
				inputFile.println("Print medium");
				inputFile.close();
			} catch(IOException e) {
			}
		}
	}
	
	public Boolean read_GronOR_arx() {
		Integer numE=0;
		Integer numC=0;
		Double eMEBF=0.0;
		Double eNOCI=0.0;
		Integer ndx=0;
		for(int i=0; i<maxMEBFs; i++ ) {
			numEnergiesGronOR[i][0]=0;
			numEnergiesGronOR[i][1]=0;
		}
		for(int i=0; i<numMEBFs; i++) {
			String fileName = projectName.trim()+mebfName[i].trim()+"_GronOR.arx";
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
								energiesGronOR[j][i][0]=eMEBF;
								numEnergiesGronOR[i][0]=j+1;
							}
						} else if(numE<=2*numC) {
							for(int j=0; j<numE; j++) {
								card=br.readLine();
								if(j>=numC) card=br.readLine();
								ndx=j;
								if(ndx>=numC) ndx=numC;
								eMEBF=Double.valueOf(card.substring(20*ndx+6,20*ndx+26)).doubleValue();
								energiesGronOR[j][i][0]=eMEBF;
								numEnergiesGronOR[i][0]=j+1;
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
								energiesGronOR[j][i][1]=eNOCI;
								numEnergiesGronOR[i][1]=j+1;
							}
							card=br.readLine();
							for(int j=0; j<numE; j++) {
								card=br.readLine();
								for(int k=0; k<numE; k++) {
									coefGronOR[j][k][i]=Double.valueOf(card.substring(20*k,20*k+20)).doubleValue();
								}
							}
						} else if(numE<=2*numC) {
							card=br.readLine();
							for(int j=0; j<7; j++) {
								eNOCI=Double.valueOf(card.substring(20*j,20*j+20)).doubleValue();
								energiesGronOR[j][i][1]=eNOCI;
								numEnergiesGronOR[i][1]=j+1;
							}
							card=br.readLine();
							for(int j=0; j<numE; j++) {
								card=br.readLine();
								for(int k=0; k<7; k++) {
									coefGronOR[j][k][i]=Double.valueOf(card.substring(20*k,20*k+20)).doubleValue();
								}
							}
							card=br.readLine();
							for(int j=7; j<numE; j++) {
								ndx=j-7;
								eNOCI=Double.valueOf(card.substring(20*ndx,20*ndx+20)).doubleValue();
								energiesGronOR[j][i][1]=eNOCI;
								numEnergiesGronOR[i][1]=j+1;
							}
							card=br.readLine();
							for(int j=0; j<numE; j++) {
								card=br.readLine();
								for(int k=7; k<numE; k++) {
									ndx=k-7;
									coefGronOR[j][k][i]=Double.valueOf(card.substring(20*ndx,20*ndx+20)).doubleValue();
								}
							}
						}
					}
		    	}
		    	br.close();
			} catch(Exception ee) {
//				System.out.println("Arx file "+fileName+" does not exist");
				return false;
			}
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
    
}
