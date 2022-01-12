package gronor;

import java.util.*;
import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;

import cformat.PrintfWriter;

public class GronOR_NWChem extends JFrame implements ActionListener, ChangeListener, WindowListener, MouseListener {

	private static final long serialVersionUID = 4L;
	
	GronOR_NWChem(String projectName){

		Boolean fileFound = false;
		Boolean converged = false;
		Double energy=0.0;
	    String card;
	    Integer numAtoms=12;
	    Double x,y,z;
	    String nwchemInput;
	    String nwchemOutput;
	    String nwchemXYZ;
		
			try {
				nwchemInput=projectName.trim()+".nw";
				nwchemOutput=projectName.trim()+".nwout";
				nwchemXYZ=projectName.trim()+".xyz";
				PrintfWriter nwFile = new PrintfWriter(new FileWriter(nwchemInput));
				nwFile.println("start fragmentA");
				nwFile.println("echo");
				nwFile.println("basis \"ao basis\" print");
				nwFile.println("* library \"def2-tzvp\"");
				nwFile.println("end");
				nwFile.println("geometry units angstrom autosym");
			    BufferedReader br = new BufferedReader(new FileReader(nwchemXYZ));
			    card=br.readLine();
			    numAtoms=Integer.parseInt(card.trim());
			    card=br.readLine();
			    for(int i=0; i<numAtoms; i++) {
			    	card=br.readLine();
			    	nwFile.println(card);
			    }
			    br.close();
				nwFile.println("end");
				nwFile.println("dft");
				nwFile.println(" xc b3lyp");
				nwFile.println("end");
				nwFile.println("driver");
				nwFile.println("end");
				nwFile.println("task dft optimize");
				nwFile.close();
				Runtime.getRuntime().exec("mv fragmentA.nw fragmentA/");
				System.out.println("In fragmentA/ run: mpirun -n 12 nwchem "+projectName+" > "+nwchemOutput);
			} catch(IOException e) {
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
