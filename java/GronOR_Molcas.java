package gronor;

import java.util.*;
import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;

import cformat.PrintfWriter;

public class GronOR_Molcas extends JFrame implements ActionListener, ChangeListener, WindowListener, MouseListener {

	private static final long serialVersionUID = 4L;
	
	GronOR_Molcas(String projectName, Integer option, Integer option2, Integer numOcc, Integer numCASe, Integer numCASo) {
		
		//option: 0=INT, 1=SCF, 2=CASSCF, 3=CASSCF+CASPT2
		//option2: 1=S0, 2=S1, 3=T1, 4=Dp, 5=Dm, 6=S2, 7=T2
		
		String[] element= new String[] {"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar"};
		Boolean fileFound = false;
		Boolean converged = false;
		Double energy=0.0;
	    String card;
	    Integer numAtoms=12;
	    Double x,y,z;
	    String previous;
	    String molcasInput = new String();
	    String molcasXYZ = new String();
	    
//		System.out.println("GronOR_Molcas called");
		
		if(option==0) {
			// Create Molcas INT input file
			try {
				molcasInput = projectName.trim()+"_INT.input";
				molcasXYZ = projectName.trim()+"_A.xyz";
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(molcasInput));
				inputFile.println("&seward");
				inputFile.println("high cholesky");
				
			    BufferedReader br = new BufferedReader(new FileReader(molcasXYZ));
			    card=br.readLine();
			    previous=card.substring(0,3).trim();
			    numAtoms=Integer.parseInt(card.trim());
			    card=br.readLine();
			    Integer count=0;
			    for(int i=0; i<numAtoms; i++) {
			    	card=br.readLine();
			    	if(!card.substring(0,3).trim().equals(previous)) {
			    		previous=card.substring(0,3).trim();
			    		if(i>0) inputFile.println("end basis set");
			    		inputFile.println("basis set");
			    		if(previous.equals("H")) inputFile.println("h.ano-s...3s3p.");
			    		if(previous.equals("C")) inputFile.println("c.ano-s...4s3p2d.");
			    		if(previous.equals("N")) inputFile.println("n.ano-s...4s3p2d.");
			    		if(previous.equals("O")) inputFile.println("o.ano-s...4s3p2d.");
			    		if(previous.equals("F")) inputFile.println("f.ano-s...4s3p2d.");
			    		if(previous.equals("Cl")) inputFile.println("cl.ano-s...4s3p2d.");
			    		if(previous.equals("Br")) inputFile.println("br.ano-s...4s3p2d.");
			    		count=0;
			    	}
			    	count++;
			    	inputFile.println(card.substring(0,2).trim()+count+" "+card.substring(3).trim()+" Angstrom");
		    		previous=card.substring(0,3).trim();
			    }
			    inputFile.println("end basis set");
			    br.close();
			    inputFile.close();
				
			} catch(IOException e) {
			}
			
		} else if(option==1) {
			// Create Molcas SCF input file
			try {
				molcasInput = projectName.trim()+"_SCF.input";
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(molcasInput));
				inputFile.println("&scf");
				inputFile.close();
			} catch(IOException e) {
			}
			
		} else {
			// Create Molcas CASSCF+CASPT2 input file
			try {
//				System.out.println("NUMOCC="+numOcc);
				Integer Inact = numOcc - numCASe/2;
				if(option2==1) {
					molcasInput = projectName.trim()+"A_S0.input";
					PrintfWriter inputFile = new PrintfWriter(new FileWriter(molcasInput));
					inputFile.println("&rasscf");
					inputFile.println("nactel");
					inputFile.println(" "+numCASe);
					inputFile.println("spin");
					inputFile.println(" 1");
					inputFile.println("inactive");
					inputFile.println(" "+Inact);
					inputFile.println("ras2");
					inputFile.println(" "+numCASo);
					inputFile.println("prwf");
					inputFile.println("  0");
					inputFile.println("prsd");
					inputFile.println(">>> COPY "+projectName.trim()+".RasOrb.1 $CurrDir/INPORB.1_1");
				    inputFile.println(">>> COPY "+projectName.trim()+".VecDet.1 $CurrDir/"+projectName.trim()+"_001.det");
					inputFile.println(">>> COPY "+projectName.trim()+".RasOrb.1 $CurrDir/INPORB.2_1");
					inputFile.println(">>> COPY "+projectName.trim()+".VecDet.1 $CurrDir/"+projectName.trim()+"_011.det");
					if(option==3) inputFile.println("&caspt2");
					inputFile.close();
				}
				if(option2==2) {
					molcasInput = projectName.trim()+"A_S1.input";
					PrintfWriter inputFile = new PrintfWriter(new FileWriter(molcasInput));
					inputFile.println("&rasscf");
					inputFile.println("nactel");
					inputFile.println(" "+numCASe);
					inputFile.println("spin");
					inputFile.println(" 1");
					inputFile.println("inactive");
					inputFile.println(" "+Inact);
					inputFile.println("ras2");
					inputFile.println(" "+numCASo);
					inputFile.println("CIRoot");
					inputFile.println(" 1 2");
					inputFile.println(" 2");
					inputFile.println("prwf");
					inputFile.println("  0");
					inputFile.println("prsd");
					inputFile.println(">>> COPY "+projectName.trim()+".RasOrb.1 $CurrDir/INPORB.1_2");
				    inputFile.println(">>> COPY "+projectName.trim()+".VecDet.2 $CurrDir/"+projectName.trim()+"_002.det");
					inputFile.println(">>> COPY "+projectName.trim()+".RasOrb.1 $CurrDir/INPORB.2_2");
					inputFile.println(">>> COPY "+projectName.trim()+".VecDet.2 $CurrDir/"+projectName.trim()+"_012.det");
					if(option==3) {
						inputFile.println("&caspt2");
						inputFile.println("Multistate= 1 2");
					}
					inputFile.close();
				}
				if(option2==3) {
					molcasInput = projectName.trim()+"A_T1.input";
					PrintfWriter inputFile = new PrintfWriter(new FileWriter(molcasInput));
					inputFile.println("&rasscf");
					inputFile.println("nactel");
					inputFile.println(" "+numCASe);
					inputFile.println("spin");
					inputFile.println(" 3");
					inputFile.println("inactive");
					inputFile.println(" "+Inact);
					inputFile.println("ras2");
					inputFile.println(" "+numCASo);
					inputFile.println("prwf");
					inputFile.println("  0");
					inputFile.println("prsd");
					inputFile.println(">>> COPY "+projectName.trim()+".RasOrb.1 $CurrDir/INPORB.1_3");
				    inputFile.println(">>> COPY "+projectName.trim()+".VecDet.1 $CurrDir/"+projectName.trim()+"_003.det");
					inputFile.println(">>> COPY "+projectName.trim()+".RasOrb.1 $CurrDir/INPORB.2_3");
					inputFile.println(">>> COPY "+projectName.trim()+".VecDet.1 $CurrDir/"+projectName.trim()+"_013.det");
					if(option==3) inputFile.println("&caspt2");
					inputFile.close();
				}
				if(option2==4) {
					molcasInput = projectName.trim()+"A_D+.input";
					PrintfWriter inputFile = new PrintfWriter(new FileWriter(molcasInput));
					inputFile.println("&rasscf");
					inputFile.println("nactel");
					inputFile.println(" "+(numCASe-1));
					inputFile.println("spin");
					inputFile.println(" 2");
					inputFile.println("inactive");
					inputFile.println(" "+Inact);
					inputFile.println("ras2");
					inputFile.println(" "+numCASo);
					inputFile.println("prwf");
					inputFile.println("  0");
					inputFile.println("prsd");
					inputFile.println(">>> COPY "+projectName.trim()+".RasOrb.1 $CurrDir/INPORB.1_4");
				    inputFile.println(">>> COPY "+projectName.trim()+".VecDet.1 $CurrDir/"+projectName.trim()+"_004.det");
					inputFile.println(">>> COPY "+projectName.trim()+".RasOrb.1 $CurrDir/INPORB.2_4");
					inputFile.println(">>> COPY "+projectName.trim()+".VecDet.1 $CurrDir/"+projectName.trim()+"_014.det");
					if(option==3) inputFile.println("&caspt2");
					inputFile.close();
				}
				if(option2==5) {
					molcasInput = projectName.trim()+"A_D-.input";
					PrintfWriter inputFile = new PrintfWriter(new FileWriter(molcasInput));
					inputFile.println("&rasscf");
					inputFile.println("nactel");
					inputFile.println(" "+(numCASe+1));
					inputFile.println("spin");
					inputFile.println(" 2");
					inputFile.println("inactive");
					inputFile.println(" "+Inact);
					inputFile.println("ras2");
					inputFile.println(" "+numCASo);
					inputFile.println("prwf");
					inputFile.println("  0");
					inputFile.println("prsd");
					inputFile.println(">>> COPY "+projectName.trim()+".RasOrb.1 $CurrDir/INPORB.1_5");
				    inputFile.println(">>> COPY "+projectName.trim()+".VecDet.1 $CurrDir/"+projectName.trim()+"_005.det");
					inputFile.println(">>> COPY "+projectName.trim()+".RasOrb.1 $CurrDir/INPORB.2_5");
					inputFile.println(">>> COPY "+projectName.trim()+".VecDet.1 $CurrDir/"+projectName.trim()+"_015.det");
					if(option==3) inputFile.println("&caspt2");
					inputFile.close();
				}
				if(option2==6) {
					molcasInput = projectName.trim()+"A_S2.input";
					PrintfWriter inputFile = new PrintfWriter(new FileWriter(molcasInput));
					inputFile.println("&rasscf");
					inputFile.println("nactel");
					inputFile.println(" "+numCASe);
					inputFile.println("spin");
					inputFile.println(" 1");
					inputFile.println("inactive");
					inputFile.println(" "+Inact);
					inputFile.println("ras2");
					inputFile.println(" "+numCASo);
					inputFile.println("CIRoot");
					inputFile.println(" 1 3");
					inputFile.println(" 3");
					inputFile.println("prwf");
					inputFile.println("  0");
					inputFile.println("prsd");
					inputFile.println(">>> COPY "+projectName.trim()+".RasOrb.1 $CurrDir/INPORB.1_6");
				    inputFile.println(">>> COPY "+projectName.trim()+".VecDet.3 $CurrDir/"+projectName.trim()+"_006.det");
					inputFile.println(">>> COPY "+projectName.trim()+".RasOrb.1 $CurrDir/INPORB.2_6");
					inputFile.println(">>> COPY "+projectName.trim()+".VecDet.3 $CurrDir/"+projectName.trim()+"_016.det");
					if(option==3) {
						inputFile.println("&caspt2");
						inputFile.println("Multistate= 1 3");
					}
					inputFile.close();
				}
				if(option2==7) {
					molcasInput = projectName.trim()+"A_T2.input";
					PrintfWriter inputFile = new PrintfWriter(new FileWriter(molcasInput));
					inputFile.println("&rasscf");
					inputFile.println("nactel");
					inputFile.println(" "+numCASe);
					inputFile.println("spin");
					inputFile.println(" 3");
					inputFile.println("inactive");
					inputFile.println(" "+Inact);
					inputFile.println("ras2");
					inputFile.println(" "+numCASo);
					inputFile.println("CIRoot");
					inputFile.println(" 1 2");
					inputFile.println(" 2");
					inputFile.println("prwf");
					inputFile.println("  0");
					inputFile.println("prsd");
					inputFile.println(">>> COPY "+projectName.trim()+".RasOrb.1 $CurrDir/INPORB.1_7");
				    inputFile.println(">>> COPY "+projectName.trim()+".VecDet.2 $CurrDir/"+projectName.trim()+"_007.det");
					inputFile.println(">>> COPY "+projectName.trim()+".RasOrb.1 $CurrDir/INPORB.2_7");
					inputFile.println(">>> COPY "+projectName.trim()+".VecDet.2 $CurrDir/"+projectName.trim()+"_017.det");
					if(option==3) {
						inputFile.println("&caspt2");
						inputFile.println("Multistate= 1 2");
					}
					inputFile.close();
				}
				
			} catch(IOException e) {
			}
			
		}
	setVisible(true);
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
