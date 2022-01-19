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

import cformat.PrintfWriter;

public class GronOR_Project extends JFrame implements ActionListener, ChangeListener, WindowListener, MouseListener{
	
	private static final long serialVersionUID = 4L;

	JMenuBar menubar = new JMenuBar();
	JMenuItem create, open, close, list, save, quit;
	String projectName, projectFile;

	JPanel parametersPanel;
	JPanel fragmentsPanel;
	JPanel fragmentsButtonsPanel;
	JPanel energiesPanel;
	JPanel mebfsPanel;
	JPanel mebfsEnergiesPanel;
	
	Integer maxFragments = 32;	// Maximum number of fragments

	/*
    numFragments = Integer.valueOf(card.substring(0,6).trim());
    newFragments = numFragments;
    dimFragments = new Integer[numFragments][11];
    namFragments = new String[numFragments];
    movFragments = new Double[numFragments][6];
    energiesDFT = new Double[numFragments];
    energiesSCF = new Double[numFragments];
    energiesCASSCF = new Double[numFragments][7];
    energiesCASPT2 = new Double[numFragments][7];
    */
	
	Integer numFragments;		// Number of fragments

	Integer[][] dimFragments = new Integer[maxFragments][11];	// Dimensions of fragments: number of atoms, states, electrons, etc.
	String[] namFragments = new String[maxFragments];			// Names of fragments used to determine xyz-formatted coordinate origin 
	Double[][] movFragments = new Double[maxFragments][6];		// Rotation and translation of coordinates with respect to original source	
	Double[] energiesDFT = new Double[maxFragments];			// DFT optimized energy of S0 state from NWChem
	Double[] energiesSCF = new Double[maxFragments];			// SCF energy of S0 state from Molcas
	Double[][] energiesCASSCF = new Double[maxFragments][7];	// CASSCF energies of all states of fragment from Molcas
	Double[][] energiesCASPT2 = new Double[maxFragments][7];	// CASPT2 energies of all states of fragment from Molcas
	
	Integer numRanks=12;		// Number of ranks in internal mpirun
	Integer numTBD=0;
	
	Integer numEnergies = 0;	// Number of energy entries

    Integer numMEBFs = 0;
    Integer newMEBFs = 0;
    Integer maxMEBFs = 32;
    Integer maxMer = 3;
    
    String[] fragmentLabels = new String[] {"ID", "SRC", "XYZ File", "Atoms", "States", "Electrons", "CASe", "CASo", "Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "Alt"};
    String[] fragmentNames = new String[] {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"};
    String[] energyNames = new String[] {"E(DFT)", "E(SCF)", "E(CASSCF)", "E(CASPT2)"};
    String[] stateLabels = new String[] {" ", "S0", "S1", "T1", "D-", "D+", "S2", "T2"};
    String[] stateNames = new String[] {"S0", "S1", "T1", "D-", "D+", "S2", "T2"};
	String[] mebfLabels = new String[] {"ID", "n-mer", "Spin", "States", "Frag", "1", "2", "3", "4", "5", "6", "7", "8", "9" };
	String[] nociLabels = new String[] {"ID", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9" };
	
    JTable fragmentsTable;
    JTable energiesTable;
    JTable mebfsTable;
    JTable mebfsEnergiesTable;
	JTable dimensionTable = new JTable();
	JTable numberTable = new JTable();
    
    Boolean numFragmentsChanged = false;

	DefaultTableModel fragmentsTableModel = new DefaultTableModel(null,fragmentLabels);
	DefaultTableModel energiesTableModel = new DefaultTableModel(null,stateLabels);
	DefaultTableModel mebfsTableModel = new DefaultTableModel(null,mebfLabels);
	DefaultTableModel mebfsEnergiesTableModel = new DefaultTableModel(null,nociLabels);

    Object[][] fragmentDefinitions = new Object[maxFragments][15];
    Object[][] stateEnergies = new Object[maxFragments][8];
    Object[][] mebfDefinitions = new Object[maxMEBFs*maxMer][14];
    Object[][] mebfEnergies = new Object[maxMEBFs][10];
    
    Integer[][] mebfIndex = new Integer[maxMEBFs*maxMer][2];
    Integer[][] mebfSpecification = new Integer[maxMEBFs][4];	 		// 0: n-mer; 1: spin; 2:number of states to include 3:fragment; 
    String[] mebfName = new String[maxMEBFs];							// name of mefb, e.g. AB for fragment A,B combination
    Integer[][][] mebfFragments = new Integer[maxMEBFs][maxMer][12];	// index to fragment states included in this mebf
    
    Integer numberStateEnergies = 0;
    Integer numberMebfDefinitions = 0;
    Integer numberMebfEnergies = 0;
    Integer numberFragmentDefinitions = 0;
    
    Integer[][] energyFragment = new Integer[maxFragments][2];

    GronOR_Fragment fragment = new GronOR_Fragment();

    Double[] RandT = new Double[6];

	
	
    JFrame dialogFrame;

    JFileChooser chooser;
    ExtensionFilter filter;
    String prjFile;
    
    
    Integer numAtoms=12, numStates=0, newFragments=0;

    
    String[] states = new String[] {" ", "S0", "S1", "T1", "D-", "D+", "S2", "T2"};
    
    String[] methods = new String[] {"SCF", "CASSCF", "CASPT2"};
    JLabel[] methodLabel = new JLabel[3];
    
    Integer[] fragments;

	
    Integer numDFT=0, numSCF=0, numCASSCF=0, numCASPT2=0;
    
    Integer numCASe=8, numCASo=8, numAlt=0, numOcc=0;
    Integer[][] alter = new Integer[12][2];
    
    
    Container container = null;

	GronOR_Project(String pName){
		
		super("GronOR Project "+pName);
	    super.setSize(1000,800);
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
	    	  if(numFragmentsChanged) System.exit(0);
	    }});

		
		if(!readProjectFile(projectFile)) initializeProject();
		
		initializeWindow();
		
	    setVisible(true);
	}

	private Boolean readProjectFile(String fileName) {

		try {
		    	BufferedReader br = new BufferedReader(new FileReader(fileName));
		    	String card;
		    	card=br.readLine();
			    numFragments = Integer.valueOf(card.substring(0,6).trim());
		    	newFragments = numFragments;
		    	numMEBFs = Integer.valueOf(card.substring(6,12).trim());
		    	newMEBFs = numMEBFs;
			    
			    for(int i=0; i<numFragments; i++) {
			    	card=br.readLine();
			    	namFragments[i]=card.trim();
			    	card=br.readLine();
			    	for(int j=0; j<6; j++) movFragments[i][j]=Double.valueOf(card.substring(j*20,j*20+20)).doubleValue();
			    	card=br.readLine();
			    	dimFragments[i][0]=Integer.valueOf(card.substring(0,6).trim());    // Index to XYZ
			    	dimFragments[i][1]=Integer.valueOf(card.substring(6,12).trim());   // Number of atoms
			    	dimFragments[i][2]=Integer.valueOf(card.substring(12,18).trim());  // Number of states
			    	dimFragments[i][3]=Integer.valueOf(card.substring(18,24).trim());  // Number of electrons
			    	dimFragments[i][4]=Integer.valueOf(card.substring(24,30).trim());  // Number of CAS electrons
			    	dimFragments[i][5]=Integer.valueOf(card.substring(30,36).trim());  // Number of CAS orbitals
			    	dimFragments[i][6]=Integer.valueOf(card.substring(36,42).trim());  // Number of DFT energy
			    	dimFragments[i][7]=Integer.valueOf(card.substring(42,48).trim());  // Number of SCF energy
			    	dimFragments[i][8]=Integer.valueOf(card.substring(48,54).trim());  // Number of CASSCF energies
			    	dimFragments[i][9]=Integer.valueOf(card.substring(54,60).trim());  // Number of CASPT2 energies
			    	dimFragments[i][10]=Integer.valueOf(card.substring(60,66).trim()); // Number of orbital alterations
			    	card=br.readLine();
			    	if(dimFragments[i][6]>0) energiesDFT[i]=Double.valueOf(card.substring(0,20)).doubleValue();
			    	if(dimFragments[i][7]>0) energiesSCF[i]=Double.valueOf(card.substring(20,40)).doubleValue();
			    	card=br.readLine();
			    	for(int j=0; j<dimFragments[i][8]; j++) energiesCASSCF[i][j]=Double.valueOf(card.substring(j*20,j*20+20)).doubleValue();
			    	card=br.readLine();
			    	for(int j=0; j<dimFragments[i][9]; j++) energiesCASPT2[i][j]=Double.valueOf(card.substring(j*20,j*20+20)).doubleValue();
			    }
			    
				for(int i=0; i<numMEBFs; i++) {
			    	card=br.readLine();
					mebfName[i]=card.trim();
			    	card=br.readLine();
					for(int j=0; j<4; j++) mebfSpecification[i][j]=Integer.valueOf(card.substring(j*6,j*6+6).trim());
					Integer nmer = mebfSpecification[i][0];
					for(int j=0; j<nmer; j++) {
				    	card=br.readLine();
						for(int k=0; k<12; k++) mebfFragments[i][j][k]=Integer.valueOf(card.substring(k*6,k*6+6).trim());
					}
				}
				
				br.close();
				return true;
			} catch(Exception ee) {
				System.out.println("Project file does not exist");
				return false;
			}
	}

	private Boolean writeProjectFile(String fileName) {
		numFragmentsChanged = numFragments!=newFragments;
		try {
			    PrintfWriter fw = new PrintfWriter(new FileWriter(fileName));
			    if(newFragments<numFragments) numFragments=newFragments;
			    fw.printf("%6d",newFragments);
			    fw.printf("%6d",numMEBFs);
			    fw.println();
				for(int i=0; i<numFragments; i++) {
					fw.println(namFragments[i].trim());
					for(int j=0; j<6; j++) fw.printf("%20.10f",movFragments[i][j]);
				    fw.println();
				    for(int j=0; j<11; j++) fw.printf("%6d",dimFragments[i][j]);
				    fw.println();
				    if(dimFragments[i][6]>0) fw.printf("%20.10f",energiesDFT[i]);
				    if(dimFragments[i][7]>0) fw.printf("%20.10f",energiesSCF[i]);
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
					for(int j=0; j<4; j++) fw.printf("%6d",mebfSpecification[i][j]); fw.println();
					Integer nmer = mebfSpecification[i][0];
					for(int j=0; j<nmer; j++) {
						for(int k=0; k<12; k++) fw.printf("%6d",mebfFragments[i][j][k]); fw.println();
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
		
		numFragments=0;
		newFragments=0;
	    numMEBFs = 0;
	    newMEBFs = 0;
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
	    	extFile = chooser.getSelectedFile().toString();
	    };
		rootName=extFile.substring(extFile.lastIndexOf("/")+1,extFile.indexOf(extension));
		return rootName;
	}

	public Integer getNumAxyz(String root) {
		String fileName=root+".xyz";
		String card;
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			card=br.readLine();
			Integer nA=Integer.valueOf(card.substring(0,6).trim());
			br.close();
			return nA;
		} catch(IOException e) {
			System.out.println("IOException in XYZ file "+fileName);
			return 0;
		}
	}

	public Integer getNumExyz(String root) {
		String fileName=root+".xyz";
		System.out.println("FILE="+fileName);
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

	private void updateFragmentList() {

		if(numFragments>numberFragmentDefinitions) {
			for(int i=numberFragmentDefinitions; i<newFragments; i++) {
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
		    	dimFragments[i][2]=5;
		    	dimFragments[i][4]=8;
		    	dimFragments[i][5]=8;
			}
		}
		
		numFragments=newFragments;
		
		for(int i=0; i<numFragments; i++) {
			fragmentDefinitions[i][0]=fragmentNames[i];
			fragmentDefinitions[i][1]=fragmentDefinitions[dimFragments[i][0]][0];
			fragmentDefinitions[i][2]=namFragments[i];
			fragmentDefinitions[i][3]=dimFragments[i][1];
			fragmentDefinitions[i][4]=dimFragments[i][2];
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
		
		fragmentsPanel.setPreferredSize(new Dimension(Short.MAX_VALUE,((numFragments)*15+55)));
		fragmentsPanel.setMinimumSize(new Dimension(Short.MAX_VALUE,((numFragments)*15+55)));
				
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
			fragmentDefinitions[i][4]=dimFragments[i][2];
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
		
		Integer maxF = -1;
		for(int i=0; i<numEnergies; i++) {
			if(dimFragments[energyFragment[i][0]][0]>maxF) maxF = dimFragments[energyFragment[i][0]][0];
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
			fragment.write_Run_Script_Fragments(i,dimFragments[i][2],numRanks);
		}
		
		Double energy = 0.0;
		Integer numAlt = 0;
		
		for(int i=0; i<=maxF; i++) {
			for(int j=0; j<6; j++) RandT[j]=movFragments[i][j];
			nameP=projectName.trim();
			nameF=namFragments[i].trim();
			nameA=" ";
			nameB=(String) fragmentDefinitions[dimFragments[i][0]][0];

			fragment.initialize(nameP, nameA, nameF, nameB, dimFragments[i][3], RandT);
			for(int j=0; j<numEnergies; j++) {
				if(i==dimFragments[energyFragment[j][0]][0]) {
					// DFT energy from NWChem optimization
					if(energyFragment[j][1]==0) {
						if(fragment.NWChem_Converged(i)) {
							energy=fragment.NWChem_DFT(i);
							energiesDFT[energyFragment[j][0]]=energy;
							dimFragments[i][6]=1;
							stateEnergies[j][1]=energy;
							energiesTableModel.setValueAt(stateEnergies[j][1],j,1);
						} else {
							fragment.write_NWChem_DFT(i);
						}
						fragment.write_NWChem_DFT(i);
					}
					// SCF energy from Molcas
					if(energyFragment[j][1]==1) {
						if(fragment.Molcas_SCF_Converged(energyFragment[j][0],dimFragments[energyFragment[j][0]][4])) {
							energy=fragment.Molcas_SCF(energyFragment[j][0],dimFragments[energyFragment[j][0]][4]);
							numAlt=fragment.Molcas_numAlt();
							dimFragments[i][10]=numAlt;
							energiesSCF[energyFragment[j][0]]=energy;
							stateEnergies[j][1]=energy;
							dimFragments[i][7]=1;
							energiesTableModel.setValueAt(stateEnergies[j][1],j,1);
						} else {
							fragment.write_Molcas_Int(energyFragment[j][0]);
							fragment.write_Molcas_SCF(energyFragment[j][0]);
						}
						fragment.write_Molcas_Int(energyFragment[j][0]);
						fragment.write_Molcas_SCF(energyFragment[j][0]);
					}
					// CASSCF energy from Molcas
					if(energyFragment[j][1]==2) {
						for(int k=0; k<dimFragments[energyFragment[j][0]][2]; k++) {
							if(fragment.Molcas_CASSCF_Converged(energyFragment[j][0],k)>=1) {
								energy=fragment.Molcas_CASSCF(energyFragment[j][0],k);
								energiesCASSCF[energyFragment[j][0]][k]=energy;
								stateEnergies[j][k]=energy;
								dimFragments[i][8]=k+1;
								energiesTableModel.setValueAt(stateEnergies[j][k],j,k+1);
							} else {
								fragment.write_Molcas_Int(energyFragment[j][0]);
								fragment.write_Molcas_SCF(energyFragment[j][0]);
								fragment.write_Molcas_CASSCF(energyFragment[j][0],k,false,dimFragments[energyFragment[j][0]][4],dimFragments[energyFragment[j][0]][5]);
							}
							fragment.write_Molcas_Int(energyFragment[j][0]);
							fragment.write_Molcas_SCF(energyFragment[j][0]);
							fragment.write_Molcas_CASSCF(energyFragment[j][0],k,false,dimFragments[energyFragment[j][0]][4],dimFragments[energyFragment[j][0]][5]);
						}
					}
					// CASPT2 energy from Molcas
					if(energyFragment[j][1]==3) {
						for(int k=0; k<dimFragments[energyFragment[j][0]][2]; k++) {
							if(fragment.Molcas_CASSCF_Converged(energyFragment[j][0],k)>=1) {
								energy=fragment.Molcas_CASPT2(energyFragment[j][0],k);
								energiesCASPT2[energyFragment[j][0]][k]=energy;
								stateEnergies[j][k]=energy;
								dimFragments[i][9]=k+1;
								energiesTableModel.setValueAt(stateEnergies[j][k],j,k+1);
							} else {
								fragment.write_Molcas_Int(energyFragment[j][0]);
								fragment.write_Molcas_SCF(energyFragment[j][0]);
								fragment.write_Molcas_CASSCF(energyFragment[j][0],k,true,dimFragments[energyFragment[j][0]][4],dimFragments[energyFragment[j][0]][5]);
							}
							fragment.write_Molcas_Int(energyFragment[j][0]);
							fragment.write_Molcas_SCF(energyFragment[j][0]);
							fragment.write_Molcas_CASSCF(energyFragment[j][0],k,true,dimFragments[energyFragment[j][0]][4],dimFragments[energyFragment[j][0]][5]);
						
						}
					}	
				}
			}
		}
		

		
		numberStateEnergies=numEnergies;
		fragmentsPanel.setPreferredSize(new Dimension(Short.MAX_VALUE,((numFragments)*15+55)));
		energiesPanel.setPreferredSize(new Dimension(Short.MAX_VALUE,((numEnergies)*20+30)));
		fragmentsPanel.setMinimumSize(new Dimension(Short.MAX_VALUE,((numFragments)*15+55)));
		energiesPanel.setMinimumSize(new Dimension(Short.MAX_VALUE,((numEnergies)*20+30)));

		fragmentsPanel.revalidate();
		fragmentsPanel.repaint();
		container.revalidate();
		container.repaint();
		container.setVisible(true);
	}

	private void updateMEBFDefinitions() {

		Integer index = 0;
		// Add spin 1 1-mers as additional MEBFs
		if(newMEBFs>numMEBFs) {
			for(int i=numMEBFs; i<newMEBFs; i++) {
				mebfSpecification[i][0]=1;		// n-mer (i.e. number of fragments, default monomer)
				mebfSpecification[i][1]=1;		// spin  (default spin 1)
				mebfSpecification[i][2]=1;		// number of states included (default one)
				mebfSpecification[i][3]=1;		// fragment index (default first)
				mebfName[i]=fragmentNames[0];	// mebf name (default name of first fragment)
				for(int j=0; j<maxMer; j++) {
					mebfFragments[i][j][0]=j;		// fragment index
					for(int k=1; k<12; k++) {
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
			for(int k=0; k<14; k++) {
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
			mebfDefinitions[index][3]=mebfSpecification[i][2];		// number of states included
			// loop of number of fragments
			for(int j=0; j<mebfSpecification[i][0]; j++) {
				mebfDefinitions[index][4]=fragmentNames[mebfFragments[i][j][0]];
				// loop over number of states defined
				for(int k=0; k<mebfSpecification[i][2]; k++) {
					mebfDefinitions[index][5+k]=stateNames[mebfFragments[i][j][k+1]];
				}
				mebfIndex[index][0]=i;
				mebfIndex[index][1]=j;
				index++;
			}
		}
		
		//Fill the Table Model
		for(int i=0; i<numRows; i++) {
			for(int j=0; j<14; j++) {
				mebfsTableModel.setValueAt(mebfDefinitions[i][j],i,j);
			}
		}

		mebfsPanel.setPreferredSize(new Dimension(Short.MAX_VALUE,((numRows)*20+35)));
		mebfsPanel.setMinimumSize(new Dimension(Short.MAX_VALUE,((numRows)*20+35)));	
		write_Molcas_MEBF_files();
		write_GronOR_NOCI();
	}

	private void updateMEBFEnergies() {

		if(numberMebfEnergies<numMEBFs) {
			for(int i=numberMebfEnergies; i<numMEBFs; i++) {
				mebfsEnergiesTableModel.addRow(mebfEnergies[i]);
			}
		}
		if(numberMebfEnergies>numMEBFs) {
			for(int i=numMEBFs; i<numberMebfEnergies; i++) {
				mebfsEnergiesTableModel.removeRow(numberMebfEnergies-1-i);
			}
		}

		numberMebfEnergies=numMEBFs;
		
		for(int i=0; i<numMEBFs; i++) {
			for(int k=0; k<10; k++) {
				mebfEnergies[i][k]=" ";
			}
			mebfEnergies[i][0]=mebfName[i];
		}

		//Fill the Table Model
		for(int i=0; i<numMEBFs; i++) {
			for(int j=0; j<10; j++) {
				mebfsEnergiesTableModel.setValueAt(mebfEnergies[i][j],i,j);
			}
		}
		
		mebfsEnergiesPanel.setPreferredSize(new Dimension(Short.MAX_VALUE,((numberMebfEnergies)*20+35)));
		mebfsEnergiesPanel.setMinimumSize(new Dimension(Short.MAX_VALUE,((numberMebfEnergies)*20+35)));
	}
	
	private void initializeWindow() {
		
		container = getContentPane();
		Box baseBox = Box.createVerticalBox();
		
		container.add(baseBox);
		 
// Parameter Panel
		
		parametersPanel = new JPanel();
		parametersPanel.setLayout(new BoxLayout(parametersPanel,BoxLayout.X_AXIS));
		parametersPanel.setPreferredSize(new Dimension(Short.MAX_VALUE,60));
		TitledBorder parametersBorder = new TitledBorder(new LineBorder(Color.black),"General Parameters");
		parametersBorder.setTitleColor(Color.black);
		parametersPanel.setBorder(parametersBorder);

		JPanel dimensionPanel = new JPanel();
		dimensionPanel.setLayout(new BoxLayout(dimensionPanel,BoxLayout.X_AXIS));
		dimensionPanel.setPreferredSize(new Dimension(100,100));
		Object[][] dimensionData = new Object[][] {
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
					value = JOptionPane.showInputDialog(jf,"Enter new number of fragments");
					if(value.length()>0) newFragments=Integer.valueOf(value);
				}
				if(dimensionTable.getSelectedRow()==1) {
					value = JOptionPane.showInputDialog(jf,"Enter new number of MEBFs");
					if(value.length()>0) newMEBFs=Integer.valueOf(value);
					updateMEBFDefinitions();
					updateMEBFEnergies();
					mebfsPanel.revalidate();
					mebfsPanel.repaint();
					mebfsEnergiesPanel.revalidate();
					mebfsEnergiesPanel.repaint();
				}
				dimensionData[0][1]=newFragments;
				dimensionData[1][1]=newMEBFs;
				dimensionPanel.repaint();
				updateFragmentList();
				updateFragmentDefinitions();
				fragmentsPanel.revalidate();
				energiesPanel.revalidate();
				fragmentsPanel.repaint();
				energiesPanel.repaint();
				updateMEBFDefinitions();
				mebfsPanel.revalidate();
				mebfsPanel.repaint();
				updateMEBFEnergies();
				mebfsEnergiesPanel.revalidate();
				mebfsEnergiesPanel.repaint();
			}
		});
		
		JPanel numberPanel = new JPanel();
		numberPanel.setLayout(new BoxLayout(numberPanel,BoxLayout.X_AXIS));
		numberPanel.setPreferredSize(new Dimension(100,100));
		Object[][] numberData = new Object[][] {
			{"Ranks", numRanks},
			{"TBD", numTBD}
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
					value = JOptionPane.showInputDialog(jf,"Enter new number of ranks");
					if(value.length()>0) numRanks=Integer.valueOf(value);
				}
				numberData[0][1]=numRanks;
				numberData[1][1]=0;
			}
		});

		JPanel buttonPanel = new JPanel();
		buttonPanel.setLayout(new BoxLayout(buttonPanel,BoxLayout.X_AXIS));
		buttonPanel.setPreferredSize(new Dimension(100,100));
		JButton updateButton = new JButton("Update");
		updateButton.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				updateFragmentList();
				updateFragmentDefinitions();
				updateMEBFDefinitions();
				updateMEBFEnergies();
			}
		});
		
		dimensionPanel.add(dimensionTable);
		numberPanel.add(numberTable);
		buttonPanel.add(updateButton);
		parametersPanel.add(dimensionPanel);
		parametersPanel.add(numberPanel);
		parametersPanel.add(buttonPanel);
		
		parametersPanel.add(Box.createHorizontalGlue());
		baseBox.add(parametersPanel);

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
						}
					}
				}
				// xyz coordinate file
				if(col==2) {
					fragmentDefinitions[row][col]=getFileRoot(".xyz");
					namFragments[row]=fragmentDefinitions[row][col].toString();
					dimFragments[row][1]=getNumAxyz(fragmentDefinitions[row][col].toString());
					dimFragments[row][3]=getNumExyz(fragmentDefinitions[row][col].toString());
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
									namFragments[j]=namFragments[row];						
								}
							}
						}
					}
					updateFragmentList();
					updateFragmentDefinitions();
					updateMEBFDefinitions();
					updateMEBFEnergies();
					fragmentsPanel.revalidate();
					fragmentsPanel.repaint();
				}
				// change number of states from { S0 S1 T1 D- D+ S2 T2 }
				if(col==4) {
					value = JOptionPane.showInputDialog(jf,"Enter number of states for fragment "+fragmentDefinitions[row][0].toString());
					if(value.length()>0) dimFragments[row][2]=Integer.valueOf(value);
				}
				// change number of CAS electrons
				if(col==6) {
					value = JOptionPane.showInputDialog(jf,"Enter number of electrons in CAS for fragment "+fragmentDefinitions[row][0].toString());
					if(value.length()>0) dimFragments[row][4]=Integer.valueOf(value);
				}
				// change number of CAS orbitals
				if(col==7) {
					value = JOptionPane.showInputDialog(jf,"Enter number of orbitals in CAS for fragment "+fragmentDefinitions[row][0].toString());
					if(value.length()>0) dimFragments[row][5]=Integer.valueOf(value);
				}
				// change translation in x
				if(col==8) {
					value = JOptionPane.showInputDialog(jf,"Enter Tx for fragment "+fragmentDefinitions[row][0].toString());
					movFragments[row][0]=Double.valueOf(value).doubleValue();
				}
				// change translation in y
				if(col==9) {
					value = JOptionPane.showInputDialog(jf,"Enter Ty for fragment "+fragmentDefinitions[row][0].toString());
					movFragments[row][1]=Double.valueOf(value).doubleValue();
				}
				// change translation in z
				if(col==10) {
					value = JOptionPane.showInputDialog(jf,"Enter Tz for fragment "+fragmentDefinitions[row][0].toString());
					movFragments[row][2]=Double.valueOf(value).doubleValue();
				}
				// change rotation in x
				if(col==11) {
					value = JOptionPane.showInputDialog(jf,"Enter Rx for fragment "+fragmentDefinitions[row][0].toString());
					movFragments[row][3]=Double.valueOf(value).doubleValue();
				}
				// change rotation in y
				if(col==12) {
					value = JOptionPane.showInputDialog(jf,"Enter Ry for fragment "+fragmentDefinitions[row][0].toString());
					movFragments[row][4]=Double.valueOf(value).doubleValue();
				}
				// change rotation in z
				if(col==13) {
					value = JOptionPane.showInputDialog(jf,"Enter Rz for fragment "+fragmentDefinitions[row][0].toString());
					movFragments[row][5]=Double.valueOf(value).doubleValue();
				}
				if(col==14) {
					for(int j=0; j<6; j++) RandT[j]=movFragments[row][j];
//					nameP=namFragments[row].trim();
					nameP=projectName.trim();
					nameF=namFragments[row].trim();
					nameA= (String) fragmentDefinitions[dimFragments[row][0]][0];
					nameB= (String) namFragments[row];
					fragment.initialize(nameP,nameA,nameF,nameB,dimFragments[row][3],RandT);
				}
				updateFragmentList();
				updateFragmentDefinitions();
				updateMEBFDefinitions();
				updateMEBFEnergies();
				fragmentsPanel.revalidate();
				fragmentsPanel.repaint();
				mebfsPanel.revalidate();
				mebfsPanel.repaint();
				mebfsEnergiesPanel.revalidate();
				mebfsEnergiesPanel.repaint();
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
		mebfsPanel.setPreferredSize(new Dimension(Short.MAX_VALUE,numberMebfDefinitions*15+150));
		mebfsPanel.setMinimumSize(new Dimension(Short.MAX_VALUE,numberMebfDefinitions*15+150));
		TitledBorder mebfsBorder = new TitledBorder(new LineBorder(Color.black),"MEBF Definitions");
		mebfsBorder.setTitleColor(Color.black);
		mebfsPanel.setBorder(mebfsBorder);
		mebfsTable = new JTable(mebfsTableModel);
		mebfsTable.getColumnModel().getColumn(0).setMaxWidth(40);
		ListSelectionModel mebfsSelectionModel = mebfsTable.getSelectionModel();
		mebfsSelectionModel.addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent e) {
				Integer row = mebfsTable.getSelectedRow();
				Integer col = mebfsTable.getSelectedColumn();
				updateMEBFDefinitions();
				updateMEBFEnergies();
				mebfsPanel.revalidate();
				mebfsPanel.repaint();
				mebfsEnergiesPanel.revalidate();
				mebfsEnergiesPanel.repaint();
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
		TableColumn fragmentColumn = mebfsTable.getColumnModel().getColumn(4);
		fragmentColumn.setCellEditor(new DefaultCellEditor(fragmentCombo));
		fragmentCombo.addMouseListener(new MouseAdapter() {
			public void mouseReleased(MouseEvent e) {
				Integer row = mebfsTable.getSelectedRow();
				Integer col = mebfsTable.getSelectedColumn();
				if(row>=0 && col==4) {
					Integer mymebf = mebfIndex[row][0];
					Integer mynmer = mebfIndex[row][1];
					mebfFragments[mymebf][mynmer][0] = fragmentCombo.getSelectedIndex();
				}
			}
		});
		
		DefaultComboBoxModel stateComboModel = new DefaultComboBoxModel(stateNames);
		JComboBox stateCombo = new JComboBox();
		stateCombo.setModel(stateComboModel);
		for(int i=0; i<9; i++) {
			TableColumn stateColumn = mebfsTable.getColumnModel().getColumn(5+i);
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
				updateMEBFDefinitions();
				updateMEBFEnergies();
				
				fragmentsPanel.revalidate();
				energiesPanel.revalidate();
				mebfsPanel.revalidate();
				fragmentsPanel.repaint();
				energiesPanel.repaint();
				mebfsPanel.repaint();
				mebfsEnergiesPanel.revalidate();
				mebfsEnergiesPanel.repaint();
				container.revalidate();
				container.repaint();
				container.setVisible(true);
			}
		});
		
		mebfsPanel.add(mebfsScroll);
		baseBox.add(mebfsPanel);

		
		mebfsEnergiesPanel = new JPanel();
		mebfsEnergiesPanel.setLayout(new BoxLayout(mebfsEnergiesPanel,BoxLayout.X_AXIS));
		mebfsEnergiesPanel.setPreferredSize(new Dimension(Short.MAX_VALUE,numberMebfEnergies*15+150));
		mebfsEnergiesPanel.setMinimumSize(new Dimension(Short.MAX_VALUE,numberMebfEnergies*15+150));
		TitledBorder mebfsEnergiesBorder = new TitledBorder(new LineBorder(Color.black),"MEBF Energies");
		mebfsEnergiesBorder.setTitleColor(Color.black);
		mebfsEnergiesPanel.setBorder(mebfsEnergiesBorder);
		mebfsEnergiesTable = new JTable(mebfsEnergiesTableModel);
		mebfsEnergiesTable.getColumnModel().getColumn(0).setMaxWidth(40);
		ListSelectionModel mebfsEnergiesSelectionModel = mebfsEnergiesTable.getSelectionModel();
		mebfsEnergiesSelectionModel.addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent e) {
				Integer row = mebfsEnergiesTable.getSelectedRow();
				Integer col = mebfsEnergiesTable.getSelectedColumn();
				mebfsEnergiesPanel.revalidate();
				mebfsEnergiesPanel.repaint();
			}
		});
		JScrollPane mebfsEnergiesScroll = new JScrollPane(mebfsEnergiesTable);
		mebfsEnergiesTable.setCellSelectionEnabled(true);
		mebfsEnergiesTable.addMouseListener(new MouseAdapter() {
			public void mouseReleased(MouseEvent e) {
					mebfEnergiesCellSelected(e);
			}
		});

		mebfsEnergiesPanel.add(mebfsEnergiesScroll);
		baseBox.add(mebfsEnergiesPanel);
		
		
		
		
		JPanel fillerPanel = new JPanel();
		fillerPanel.setLayout(new BoxLayout(fillerPanel,BoxLayout.Y_AXIS));
		TitledBorder fillerBorder = new TitledBorder(new LineBorder(Color.black),"Filler Panel");
		fillerBorder.setTitleColor(Color.black);
		fillerPanel.setBorder(fillerBorder);
        fillerPanel.add(new Box.Filler(new Dimension(10,10), new Dimension(Short.MAX_VALUE,300), new Dimension(Short.MAX_VALUE,Short.MAX_VALUE)));
		baseBox.add(fillerPanel);
		baseBox.add(Box.createVerticalGlue());

		updateFragmentList();
		updateFragmentDefinitions();
		updateMEBFDefinitions();
		updateMEBFEnergies();
		fragmentsPanel.revalidate();
		energiesPanel.revalidate();
		mebfsPanel.revalidate();
		fragmentsPanel.repaint();
		energiesPanel.repaint();
		mebfsPanel.repaint();
		mebfsEnergiesPanel.revalidate();
		mebfsEnergiesPanel.repaint();
		container.revalidate();
		container.repaint();
		container.setVisible(true);
	}
		
	private void energyCellSelected(MouseEvent e) {
	}
	
	private void mebfCellSelected(MouseEvent e) {
		Integer row = mebfsTable.getSelectedRow();
		Integer col = mebfsTable.getSelectedColumn();
		JFrame jf = new JFrame();
		String value;
		Integer mebf = mebfIndex[row][0];
		Integer nmer = mebfSpecification[mebf][0];
		Integer spin = mebfSpecification[mebf][1];
		Integer stat = mebfSpecification[mebf][2];
		if(col==1) {
			value = JOptionPane.showInputDialog(jf,"Enter number of fragments for MEBF "+mebfName[mebf].trim());
			if(value.length()>0) mebfSpecification[mebf][0]=Integer.valueOf(value);
			if(mebfSpecification[mebf][0]>maxMer) mebfSpecification[mebf][0]=maxMer;
			if(mebfSpecification[mebf][0]>nmer) {
				for(int i=nmer; i<mebfSpecification[mebf][0]; i++) {
					for(int k=0; k<12; k++) {
						mebfFragments[mebf][i][k]=0;
					}
					for(int k=0; k<mebfSpecification[mebf][3]+1; k++) {
						mebfFragments[mebf][i][k]=mebfFragments[mebf][nmer-1][k];
					}
					mebfFragments[mebf][i][0]=mebfFragments[mebf][i-1][0]+1;
					if(mebfFragments[mebf][i][0]>maxMer) mebfFragments[mebf][i][0]=mebfFragments[mebf][i][0]-maxMer;
				}
			}
			updateMEBFDefinitions();
			updateMEBFEnergies();
			mebfsPanel.revalidate();
			mebfsPanel.repaint();
			mebfsEnergiesPanel.revalidate();
			mebfsEnergiesPanel.repaint();
			container.revalidate();
			container.repaint();
			container.setVisible(true);
		}
		
		if(col==2) {
			value = JOptionPane.showInputDialog(jf,"Enter spin for MEBF "+mebfName[mebf].trim());
			if(value.length()>0) mebfSpecification[mebf][1]=Integer.valueOf(value);
			if(mebfSpecification[mebf][1]!=spin) {
				spin=mebfSpecification[mebf][1];
				if(mebfSpecification[mebf][0]==1 && mebfSpecification[mebf][1]==1) {
					for(int i=0; i<stat; i++) mebfFragments[mebf][0][i+1]=i;
				}
				if(mebfSpecification[mebf][0]==2 && mebfSpecification[mebf][1]==1) {
					mebfFragments[mebf][0][1]=0; mebfFragments[mebf][1][1]=0;
					mebfFragments[mebf][0][2]=1; mebfFragments[mebf][1][2]=0;
					mebfFragments[mebf][0][3]=0; mebfFragments[mebf][1][3]=1;
					mebfFragments[mebf][0][4]=2; mebfFragments[mebf][1][4]=2;
					mebfFragments[mebf][0][5]=3; mebfFragments[mebf][1][5]=4;
					mebfFragments[mebf][0][6]=4; mebfFragments[mebf][1][6]=3;
					mebfFragments[mebf][0][7]=1; mebfFragments[mebf][1][7]=1;
				}
				if(mebfSpecification[mebf][0]==2 && mebfSpecification[mebf][1]==1) {
					mebfFragments[mebf][0][1]=0; mebfFragments[mebf][1][1]=0;
					mebfFragments[mebf][0][2]=1; mebfFragments[mebf][1][2]=0;
					mebfFragments[mebf][0][3]=0; mebfFragments[mebf][1][3]=1;
					mebfFragments[mebf][0][4]=2; mebfFragments[mebf][1][4]=2;
					mebfFragments[mebf][0][5]=3; mebfFragments[mebf][1][5]=4;
					mebfFragments[mebf][0][6]=4; mebfFragments[mebf][1][6]=3;
					mebfFragments[mebf][0][7]=1; mebfFragments[mebf][1][7]=1;
				}
			}
			updateMEBFDefinitions();
			updateMEBFEnergies();
			mebfsPanel.revalidate();
			mebfsPanel.repaint();
			mebfsEnergiesPanel.revalidate();
			mebfsEnergiesPanel.repaint();
			container.revalidate();
			container.repaint();
			container.setVisible(true);
		}
		
		if(col==3) {
			value = JOptionPane.showInputDialog(jf,"Enter number of states for MEBF "+mebfName[mebf].trim());
			if(value.length()>0) mebfSpecification[mebf][2]=Integer.valueOf(value);
			if(mebfSpecification[mebf][2]!=stat) {
				stat=mebfSpecification[mebf][2];
				if(mebfSpecification[mebf][0]==1 && mebfSpecification[mebf][1]==1) {
					for(int i=0; i<stat; i++) mebfFragments[mebf][0][i+1]=i;
				}
				if(mebfSpecification[mebf][0]==2 && mebfSpecification[mebf][1]==1) {
					mebfFragments[mebf][0][1]=0; mebfFragments[mebf][1][1]=0;
					mebfFragments[mebf][0][2]=1; mebfFragments[mebf][1][2]=0;
					mebfFragments[mebf][0][3]=0; mebfFragments[mebf][1][3]=1;
					mebfFragments[mebf][0][4]=2; mebfFragments[mebf][1][4]=2;
					mebfFragments[mebf][0][5]=3; mebfFragments[mebf][1][5]=4;
					mebfFragments[mebf][0][6]=4; mebfFragments[mebf][1][6]=3;
					mebfFragments[mebf][0][7]=1; mebfFragments[mebf][1][7]=1;
				}
			}
			updateMEBFDefinitions();
			updateMEBFEnergies();
			mebfsPanel.revalidate();
			mebfsPanel.repaint();
			mebfsEnergiesPanel.revalidate();
			mebfsEnergiesPanel.repaint();
			container.revalidate();
			container.repaint();
			container.setVisible(true);
		}
	}

	private void mebfEnergiesCellSelected(MouseEvent e) {
	}
	
	private void write_Molcas_MEBF_files() {
		for(int i=0; i<numMEBFs; i++) {
			Integer nfrags = mebfSpecification[i][0];
			String pName = " ";
			String[] frags = new String[nfrags];
			Integer[] nums = new Integer[nfrags];
			Double[][] randt = new Double[nfrags][6];
			String fileName = projectName.trim()+mebfName[i].trim();
			if(nfrags==1) fileName=fileName.trim()+mebfName[i].trim();
			for(int j=0; j<nfrags; j++) {
				frags[j]=fragmentNames[mebfFragments[i][j][0]];
				nums[j]=dimFragments[j][2];
				pName = projectName.trim();
				for(int k=0; k<6; k++) randt[j][k]=movFragments[mebfFragments[i][j][0]][k];
			}
			if(fragment.write_MEBF_XYZ(fileName, pName, nfrags, frags)) {
				fragment.write_Molcas_MEBF_One(fileName, pName, nfrags, frags);	
				fragment.write_Molcas_MEBF_CB(fileName, nfrags, nums);	
				fragment.write_Molcas_MEBF_Two(fileName, pName, nfrags, frags);
				fragment.write_Run_Script_MEBFs(fileName,numRanks);
			}
		}
	}

	public void write_GronOR_NOCI() {
		Integer numME =0;
		Integer nmer = 0;
		for(int i=0; i<numMEBFs; i++) {
			String fileName = projectName.trim()+mebfName[i].trim()+"_GronOR.inp";
			numME=mebfSpecification[i][2];
			nmer=mebfSpecification[i][0];
			try {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("MEBFs "+projectName.trim()+" "+numME);
				for(int j=0; j<nmer; j++) {
					inputFile.print(fragmentNames[mebfFragments[i][j][0]]+" ");
					for(int k=0; k<numME; k++) inputFile.print(" "+stateNames[mebfFragments[i][j][k+1]]);
					inputFile.println();
				}
				inputFile.println("Threshold 1.0e-5");
				inputFile.println("Print medium");
				inputFile.close();
			} catch(IOException e) {
			}
		}
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
