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

public class gronor_Analysis extends JFrame implements ActionListener, ChangeListener, WindowListener, MouseListener, ItemListener{

	private static final long serialVersionUID = 4L;

	String projectName;
	String projectFile;
	
	JMenuBar menubar = new JMenuBar();
	JMenuItem create, open, close, list, quit;

	JPanel nociResultsPanel;
	JPanel nociEnergiesPanel;
	JPanel nociPlotPanel;

    Plot energyDiagram = new Plot();
    
    Container container = null;

    Integer noci0=0;
    Integer noci1=0;
    Integer numMEBFs = 0;
    Integer maxMEBFs = 60;
    Integer maxMEBFstates=60;
    Integer maxMer = 10;
    String[] mebfName = new String[maxMEBFs];				

    Integer[][] mebfSpecification = new Integer[maxMEBFs][5];	 	        	// 0: n-mer; 1: spin; 2: charge; 3:number of states to include 4:fragment; 

    JTable nociEnergiesTable;
	DefaultTableModel nociEnergiesTableModel = new DefaultTableModel();

    Object[][] nociEnergies = new Object[maxMEBFstates][maxMEBFs];

	Integer[][] numEnergiesGronOR = new Integer[maxMEBFs][2];		// Number of GronORenergy entries
	Double[][][] energiesGronOR = new Double[35][maxMEBFs][2];
	Double[][][] energiesRelGronOR = new Double[35][maxMEBFs][2];
	Double[][][] coefGronOR = new Double[35][35][maxMEBFs];

	String account="chm154";
	String jobName="GronOR";
	String timeLimit="00:30:00";
	
	Integer numSets=0, newSets, numFragments, newFragments;
	Integer newMEBFs, numRanks, numThreads, basisSet, contract, cholesky, expansion;
	Double thresh_MO, thresh_CI, ipea, fieldX, fieldY, fieldZ;
	
	Integer maxSets = 6;
	Integer maxFragments = 32;
	
	Integer[] lenStateList = new Integer[maxSets];
	Integer[] spinStateList = new Integer[maxSets];
	Integer[][] ndxStateList = new Integer[maxSets][16];
	
	Integer[][] dimFragments = new Integer[maxFragments][12];
	String[] namFragments = new String[maxFragments];
	Double[][] movFragments = new Double[maxFragments][6];
	
	Object[][] fragmentDefinitions = new Object[maxFragments][16];
	
	Double[] energiesDFT = new Double[maxFragments];
	Double[] energiesSCF = new Double[maxFragments];
	Double[][] energiesCASSCF = new Double[maxFragments][12];
	Double[][] energiesCASPT2 = new Double[maxFragments][12];
	
	Integer[][][] mebfFragments = new Integer[maxMEBFs][maxMer][maxMEBFstates];
	
	gronor_Analysis(String p){

		super("GronOR Analysis "+p.trim());
	    super.setSize(1200,900);
	    super.setLocation(100,100);
	    super.setBackground(Color.lightGray);
	    super.setForeground(Color.blue);
	    super.addWindowListener(this);
	    super.setFont(new Font("Serif",Font.BOLD,10));

	    this.setJMenuBar(menubar);

		projectName=p;
		projectFile=projectName+".prj";
		
	    JMenu select = new JMenu("File");
	    menubar.add(select);

	    select.add(close = new JMenuItem("Close Analysis Window"));
	    close.addActionListener(new ActionListener(){
	      public void actionPerformed(ActionEvent e){ 
		setVisible(false); }});
	    
	    select.add(quit = new JMenuItem("Quit GronOR GUI"));
	    quit.addActionListener(new ActionListener(){
	      public void actionPerformed(ActionEvent e){ 
		setVisible(false); System.exit(0); }});

		if(readProjectFile(projectFile)) {
			initializeWindow();
			updateNOCIEnergies();
			nociResultsPanel.setPreferredSize(new Dimension(Short.MAX_VALUE,512));
			nociResultsPanel.setMinimumSize(new Dimension(Short.MAX_VALUE,512));
			nociResultsPanel.setMaximumSize(new Dimension(Short.MAX_VALUE,noci1*15+65));
			nociEnergiesPanel.setPreferredSize(new Dimension(noci0*80+40,noci1*15+50));
			nociEnergiesPanel.setMinimumSize(new Dimension(noci0*80+40,noci1*15+50));
			nociEnergiesPanel.setMaximumSize(new Dimension(noci0*80+40,noci1*15+50));
			nociPlotPanel.setPreferredSize(new Dimension(300,noci1*15+50));
			nociPlotPanel.setMinimumSize(new Dimension(300,noci1*15+50));
			nociPlotPanel.setMaximumSize(new Dimension(300,noci1*15+50));
			nociEnergiesPanel.revalidate();
			nociPlotPanel.revalidate();
			nociResultsPanel.revalidate();
			nociEnergiesPanel.repaint();
			nociPlotPanel.repaint();
			nociResultsPanel.repaint();
			container.revalidate();
			container.repaint();
			container.setVisible(true);
		} else {
			
		}
		
	    setVisible(true);
	    
	    super.setVisible(true);
	    
	}

	
	private void initializeWindow() {
		
		container = getContentPane();
		Box baseBox = Box.createVerticalBox();
		
		container.add(baseBox);
		
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
		
		container.revalidate();
		container.repaint();
		container.setVisible(true);
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
		    	fieldX    = Double.valueOf(card.substring(100,120)).doubleValue();
		    	fieldY    = Double.valueOf(card.substring(120,140)).doubleValue();
		    	fieldZ    = Double.valueOf(card.substring(140,160)).doubleValue();
			    
		    	for(int i=0; i<numSets; i++) {
		    		card=br.readLine();
		    		lenStateList[i]=Integer.valueOf(card.substring(0,6).trim());
		    		spinStateList[i]=1;
		    		for(int j=0; j<lenStateList[i]; j++) {
		    			ndxStateList[i][j]=Integer.valueOf(card.substring((j+1)*6,(j+2)*6).trim());
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
			    	dimFragments[i][11]=Integer.valueOf(card.substring(66,72).trim()); // Charge
			    	fragmentDefinitions[i][1]=card.substring(77,78).trim(); // Source fragment number
			    	card=br.readLine();
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
	
	private void nociEnergiesCellSelected(MouseEvent e) {
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
	