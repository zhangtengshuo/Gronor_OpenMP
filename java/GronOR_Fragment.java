package gronor;

import java.util.*;
import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;

import cformat.PrintfWriter;

public class GronOR_Fragment {

	private static final long serialVersionUID = 4L;

	String projectName;
	
	Integer maxAtoms = 100;
	
	String framentName;
	String projectRoot;
	
	String[] fragmentNames = new String[] {"A", "B", "C", "D", "E", "F", "G", "H"};
    String[] stateNames = new String[] {"S0","S1","S2","D0","D1","T1","T2","S+","D+","T+","S-","D-","T-"};
	
	Integer fragmentID = 0, numAtoms=0, numOcc=0;
	String fragmentName = " ";
	String fragmentRoot = " ";
	String fragmentAXYZ, fragmentXYZ, nwchemOutput;

	Double[][] coordinates = new Double[maxAtoms][3];
	String[] atomLabel = new String[maxAtoms];

	Double[] v = new Double[] {0.0,0.0,0.0};
	Double[] w = new Double[] {0.0,0.0,0.0};
	Double[] x = new Double[] {0.0,0.0,0.0};
	Double[] y = new Double[] {0.0,0.0,0.0};
	Double g = 0.0;
	
	Double RandT[] = new Double[6];
	
	Integer numAlt;
	Integer[][] alter = new Integer[12][2];

	GronOR_Fragment(){
	}

	public void initialize(String nameP,String nameA, String nameF,String nameB, Integer numE, Double[] rt) {
		// Initialize PA from source FB
		Boolean writeFile = false; 
		projectRoot=nameP;
		fragmentRoot=nameF;
		
		if(!projectRoot.trim().equals(fragmentRoot.trim())) writeFile=true;
		if(!nameB.trim().equals(nameA.trim())) writeFile=true;
		
		fragmentName=nameF.trim()+nameB.trim();
		numOcc = numE / 2;
		for(int i=0; i<6; i++) RandT[i]=rt[i];
		
		if(read_XYZ()) {
			fragmentName=nameP.trim()+nameA.trim();
			if(writeFile) {
				rotate_AND_translate(RandT);
				write_XYZ();
			}
		} else {
			if(NWChem_Converged(0)) {
				if(!nameB.trim().equals(nameA.trim())) {
					rotate_AND_translate(RandT);
					if(!write_XYZ()) System.exit(0);
				}
			}
		}
	}

	public void initialize2(String nameF, String nameA, Integer numE) {
		fragmentRoot=nameF;
		fragmentName=nameF.trim()+nameA.trim();
		numOcc = numE / 2;
		if(read_XYZ()) {
		} else {
			System.out.println("XYZ file "+fragmentName+".xyz not found"); System.exit(0);
		}
	}

	public void write_Run_Script_Fragments(Integer frag, Integer stateset, Integer[] lenStateList, Integer[][] ndxStateList, Integer ranks) {
		String rootName = projectRoot.trim()+fragmentNames[frag].trim();
		String fileName = rootName.trim()+".run";
		String fullName;
		Integer numStates = lenStateList[stateset];
		try {
			PrintfWriter runFile = new PrintfWriter(new FileWriter(fileName));
			runFile.println("#!/usr/bin/tcsh");
			runFile.println("setenv MOLCAS_NPROCS "+ranks);
			runFile.println("setenv MOLCAS_MEM 1024");
			fullName = rootName.trim()+"_INT";
			runFile.println("cp "+fullName.trim()+".input "+rootName.trim()+".input; "+"pymolcas "+rootName.trim()+".input > "+fullName.trim()+".output");
			fullName = rootName.trim()+"_SCF";
			runFile.println("cp "+fullName.trim()+".input "+rootName.trim()+".input; "+"pymolcas "+rootName.trim()+".input > "+fullName.trim()+".output");
			for(int i=0; i<numStates; i++) {
				Integer index = ndxStateList[stateset][i];
				fullName = rootName.trim()+"_"+stateNames[index].trim();
				runFile.println("cp "+fullName.trim()+".input "+rootName.trim()+".input; "+"pymolcas "+rootName.trim()+".input > "+fullName.trim()+".output");
			}
			runFile.close();
		} catch(IOException ei) {
		}
	}

	public void write_Run_Script_MEBFs(String p, String pn, Integer nfrags, String[] frags, Integer[] states, Integer[][] ndxStateList, Integer ranks) {
		String fileName = p+"_Molcas.run";
		String fullName;
		String fullName2;
		try {
			PrintfWriter runFile = new PrintfWriter(new FileWriter(fileName));
			runFile.println("#!/usr/bin/tcsh");
			runFile.println("setenv MOLCAS_NPROCS "+ranks);
			runFile.println("setenv MOLCAS_MEM 1024");
			fullName = p.trim()+"_MEBFONE";
			runFile.println("cp "+fullName.trim()+".input "+p.trim()+".input; "+"pymolcas "+p.trim()+".input > "+fullName.trim()+".output");
			Integer index=0;
			Integer stateIndex = 0;
			for(int i=0; i<nfrags; i++) {
				runFile.println("cp "+pn.trim()+frags[i].trim()+".RUNFIL RUNFIL"+(i+1));
				runFile.println("cp "+pn.trim()+frags[i].trim()+".ONEINT ONEINT"+(i+1));
				for(int j=0; j<states[i]; j++) {
					index++;
					stateIndex=ndxStateList[i][j];
					runFile.println("cp "+pn.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".INPORB INPORB."+(i+1)+"_"+(j+1));
					if(index<10) {
						runFile.println("cp "+pn.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".det "+p.trim()+"_00"+index+".det");
					} else {
						runFile.println("cp "+pn.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".det "+p.trim()+"_0"+index+".det");
					}
				}
			}
			runFile.println("cp "+p.trim()+".RUNFILE RUNFILE");
			fullName = p.trim()+"_MEBFCB";
			runFile.println("cp "+fullName.trim()+".input "+p.trim()+".input; "+"common_basis < "+p.trim()+".input > "+fullName.trim()+".output");
			index=0;
			for(int i=0; i<nfrags; i++) {
				runFile.println("rm RUNFIL"+(i+1));
//				runFile.println("rm ONEINT"+(i+1));
				for(int j=0; j<states[i]; j++) {
					index++;
					stateIndex=ndxStateList[i][j];
					if(index<10) {
						runFile.println("mv "+p.trim()+"_00"+index+".vec "+p.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".vec");
						runFile.println("mv "+p.trim()+"_00"+index+".det "+p.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".det");
					} else {
						runFile.println("mv "+p.trim()+"_0"+index+".vec "+p.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".vec");
						runFile.println("mv "+p.trim()+"_0"+index+".det "+p.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".det");
					}
					runFile.println("rm INPORB."+(i+1)+"_"+(j+1));
				}	
			}
			fullName = p.trim()+"_MEBFCB";
			runFile.println("setenv DELETED ` grep \"Deleted orbitals in MOTRA\" "+fullName.trim()+".output | cut -b42- `");
			runFile.println("touch TRAINT");
			runFile.println("setenv MOLCAS_MEM 4096");
			fullName = p.trim()+"_MEBFTWO";
			runFile.println("cp "+fullName.trim()+".input "+p.trim()+".input; "+"pymolcas "+p.trim()+".input > "+fullName.trim()+".output");
			runFile.println("setenv OMP_NUM_THREADS 12");
			runFile.println("cp "+p.trim()+".RUNFILE RUNFILE");
			runFile.println("cp "+p.trim()+"_CHMOT1 _CHMOT1");
			fullName = p.trim()+"_MEBFRT";
			runFile.println("rdcho $MOLCAS_NPROCS > "+fullName.trim()+".output");
			runFile.println("rm _CHMOT1");
			runFile.println("cp "+p.trim()+".RUNFILE RUNFILE");
			runFile.println("cp "+p.trim()+".CHORST CHORST");
			runFile.println("cp "+p.trim()+".ONEINT ONEINT");
			runFile.println("cp "+p.trim()+".TRAONE TRAONE");
			runFile.println("cp "+p.trim()+".CHOMAP CHOMAP");
			runFile.println("cp "+p.trim()+".CHRED  CHRED");
			runFile.println("cp "+p.trim()+".CHVEC1 CHVEC1");
			runFile.println("cp "+p.trim()+".COMMONORB COMMONORB");
			fullName = p.trim()+"_MEBFCB";
			fullName2 = p.trim()+"_MEBFRT";
			runFile.println("rdtraint < "+fullName.trim()+".input > "+fullName2+".output");
			runFile.println("mv COMMONORB "+p.trim()+".COMMONORB");
			runFile.println("rm RUNFILE");
			runFile.println("rm CHORST");
			runFile.println("rm ONEINT");
			runFile.println("rm TRAONE");
			runFile.println("rm TRAINT");
			runFile.println("rm CHOMAP");
			runFile.println("rm CHRED");
			runFile.println("rm CHVEC1");
			runFile.println("unsetenv DELETED");
			runFile.close();
		} catch(IOException ei) {
		}
		fileName = p+"_GronOR.run";
		try {
			PrintfWriter runFile = new PrintfWriter(new FileWriter(fileName));
			runFile.println("#!/usr/bin/tcsh");
			fullName= p+"_GronOR";
			runFile.println("mpirun -n "+ranks+" gronor "+fullName);
			runFile.close();
		} catch(IOException ei) {
		}
	}
	
	public Boolean read_XYZ() {
		String fileName = fragmentName+".xyz";
		String card;
		StringTokenizer st;
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			card=br.readLine();
			st = new StringTokenizer(card," ");
			numAtoms=Integer.valueOf(st.nextToken());
			card=br.readLine();
			for(int i=0; i<numAtoms; i++) {
				card=br.readLine();
				st = new StringTokenizer(card," ");
				atomLabel[i]=st.nextToken();
				coordinates[i][0]=Double.valueOf(st.nextToken()).doubleValue();
				coordinates[i][1]=Double.valueOf(st.nextToken()).doubleValue();
				coordinates[i][2]=Double.valueOf(st.nextToken()).doubleValue();
			}
			br.close();
			return true;
		} catch(IOException ef) {
			return false;
		}
	}
	
	public Boolean write_XYZ() {
		String fileName = fragmentName+".xyz";
		try {
			PrintfWriter xyzFile = new PrintfWriter(new FileWriter(fileName));
			xyzFile.printf("%6d",numAtoms);
			xyzFile.println();
			xyzFile.println("Coordinates in Angstrom");
			for(int i=0; i<numAtoms; i++) {
				xyzFile.print(atomLabel[i]); 
				xyzFile.printf("%12.8f",coordinates[i][0]);
				xyzFile.printf("%12.8f",coordinates[i][1]);
				xyzFile.printf("%12.8f",coordinates[i][2]);
				xyzFile.println();
			}
			xyzFile.close();
			return true;
		} catch(IOException ei) {
			return false;
		}
	}

	public Boolean NWChem_Converged(Integer frag) {
		String fileName = projectRoot+fragmentNames[frag]+"_DFT.nwout";
		String card;
		Boolean converged = false;
		Double energy = 0.0;
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			while((card=br.readLine()) != null) {
				if(card.contains("Optimization converged")) converged=true;
				if(converged && card.startsWith("@")) energy=Double.valueOf(card.substring(7,22).trim()).doubleValue();
				if(converged && card.startsWith(" Output coordinates")) {
					card=br.readLine(); card=br.readLine(); card=br.readLine();
					try {
						
						for(int i=0; i<numAtoms; i++) {
							card=br.readLine();
							coordinates[i][0]=Double.valueOf(card.substring(34,47)).doubleValue();
							coordinates[i][1]=Double.valueOf(card.substring(49,62)).doubleValue();
							coordinates[i][2]=Double.valueOf(card.substring(64,77)).doubleValue();
						}
					} catch(IOException ei) {
					}
				}
			}
			br.close();
		} catch(IOException ef) {
			converged=false;
		}
		return converged;
	}
	
	public Double NWChem_DFT(Integer frag) {
		String fileName = projectRoot+fragmentNames[frag]+"_DFT.nwout";
		String card;
		Boolean converged = false;
		Double energy = 0.0;
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			while((card=br.readLine()) != null) {
				if(card.contains("Optimization converged")) converged=true;
				if(converged && card.startsWith("@")) energy=Double.valueOf(card.substring(7,22).trim()).doubleValue();
				if(converged && card.startsWith(" Output coordinates")) {
					card=br.readLine(); card=br.readLine(); card=br.readLine();
					try {
						
						for(int i=0; i<numAtoms; i++) {
							card=br.readLine();
							coordinates[i][0]=Double.valueOf(card.substring(34,47)).doubleValue();
							coordinates[i][1]=Double.valueOf(card.substring(49,62)).doubleValue();
							coordinates[i][2]=Double.valueOf(card.substring(64,77)).doubleValue();
						}
					} catch(IOException ei) {
					}
				}
			}
			br.close();
		} catch(IOException ef) {
			converged=false;
		}
		return energy;
	}
	
	public void write_NWChem_DFT(Integer frag, Integer ranks) {
		String fileName = projectRoot+fragmentNames[frag]+"_DFT.nw";
		try {
			PrintfWriter nwFile = new PrintfWriter(new FileWriter(fileName));
			nwFile.println("start "+projectRoot+fragmentNames[frag]);
			nwFile.println("echo");
			nwFile.println("basis \"ao basis\" print");
			nwFile.println("* library \"def2-tzvp\"");
			nwFile.println("end");
			nwFile.println("geometry units angstrom autosym");
		    for(int i=0; i<numAtoms; i++) {
		    	nwFile.print(atomLabel[i]+" "); 
				nwFile.printf("%12.8f",coordinates[i][0]);
				nwFile.printf("%12.8f",coordinates[i][1]);
				nwFile.printf("%12.8f",coordinates[i][2]);
				nwFile.println();
		    }
		    nwFile.println("end");
			nwFile.println("dft");
			nwFile.println(" xc b3lyp");
			nwFile.println("end");
			nwFile.println("driver");
			nwFile.println("end");
			nwFile.println("task dft optimize");
			nwFile.close();
		} catch(IOException e) {
		}
		fileName = projectRoot+fragmentNames[frag]+"_DFT.run";
		try {
			PrintfWriter runFile = new PrintfWriter(new FileWriter(fileName));
		    runFile.println("mpirun -n "+ranks+" nwchem "+projectRoot+fragmentNames[frag]+"_DFT > "+projectRoot+fragmentNames[frag]+"_DFT.nwout");
			runFile.close();
		} catch(IOException e) {
		}		
	}
	
	public Boolean write_Molcas_Int(Integer frag) {
		String fileName = projectRoot+fragmentNames[frag]+"_INT.input";
		String rootName=projectRoot.trim()+fragmentNames[frag].trim();
		String previous;
		try {
			PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
			inputFile.println("&seward");
			inputFile.println("high cholesky");			
		    previous=" ";
		    Integer count=0;
		    for(int i=0; i<numAtoms; i++) {
		    	if(!atomLabel[i].trim().equals(previous)) {
		    		previous=atomLabel[i].trim();
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
		    	inputFile.print(atomLabel[i].trim()+count+" "); 
				inputFile.printf("%12.8f",coordinates[i][0]);
				inputFile.printf("%12.8f",coordinates[i][1]);
				inputFile.printf("%12.8f",coordinates[i][2]);
				inputFile.println(" Angstrom");
				previous=atomLabel[i].trim();
		    }
		    inputFile.println("end basis set");
			inputFile.println(">>> COPY "+rootName.trim()+".OneInt $CurrDir/"+rootName.trim()+".ONEINT");
			inputFile.println(">>> COPY "+rootName.trim()+".RunFile $CurrDir/"+rootName.trim()+".RUNFIL");
		    inputFile.close();
		    return true;
		} catch(IOException e) {
			return false;
		}
	}
	
	public Boolean write_Molcas_SCF(Integer frag) {
		String fileName = projectRoot+fragmentNames[frag]+"_SCF.input";
		try {
			PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
			inputFile.println("&scf");
			inputFile.close();
			return true;
		} catch(IOException e) {
			return false;
		}
	}
	
	public Integer Molcas_numAlt() {
		return numAlt;
	}
	
	public Boolean write_Molcas_CASSCF(Integer frag, Integer state, Boolean withCASPT2, Integer numCASe, Integer numCASo) {
		String fileName = projectRoot.trim()+fragmentNames[frag].trim()+"_"+stateNames[state].trim()+".input";
		String rootName=projectRoot.trim()+fragmentNames[frag].trim();
		String ext = "_"+stateNames[state];
		System.out.println("IN WRITE_MOLCAS for "+fragmentNames[frag].trim()+" "+stateNames[state].trim()+" : "+fileName+" "+state);
		try {
			Integer Inact = numOcc - numCASe/2;
			if(stateNames[state].trim().equals("S0")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
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
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				if(withCASPT2) inputFile.println("&caspt2");
				inputFile.close();
				return true;
			} else if(stateNames[state].trim().equals("S1")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
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
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.2 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("Multistate= 1 2");
				}
				inputFile.close();
				return true;
			} else if(stateNames[state].trim().equals("S2")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
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
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.3 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("Multistate= 1 2");
				}
				inputFile.close();
				return true;
			} else if(stateNames[state].trim().equals("D0")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				inputFile.println("nactel");
				inputFile.println(" "+numCASe);
				inputFile.println("spin");
				inputFile.println(" 2");
				inputFile.println("inactive");
				inputFile.println(" "+Inact);
				inputFile.println("ras2");
				inputFile.println(" "+numCASo);
				inputFile.println("prwf");
				inputFile.println("  0");
				inputFile.println("prsd");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.2 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("Multistate= 1 2");
				}
				inputFile.close();
				return true;
			} else if(stateNames[state].trim().equals("D1")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				inputFile.println("nactel");
				inputFile.println(" "+numCASe);
				inputFile.println("spin");
				inputFile.println(" 2");
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
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.2 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("Multistate= 1 2");
				}
				inputFile.close();
				return true;
			} else if(stateNames[state].trim().equals("T1")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
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
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				if(withCASPT2) inputFile.println("&caspt2");
				inputFile.close();
				return true;
			} else if(stateNames[state].trim().equals("T2")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
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
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.2 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("Multistate= 1 2");
				}
				inputFile.close();
				return true;
			} else if(stateNames[state].trim().equals("S+")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				inputFile.println("nactel");
				inputFile.println(" "+(numCASe-1));
				inputFile.println("spin");
				inputFile.println(" 1");
				inputFile.println("inactive");
				inputFile.println(" "+Inact);
				inputFile.println("ras2");
				inputFile.println(" "+numCASo);
				inputFile.println("prwf");
				inputFile.println("  0");
				inputFile.println("prsd");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				if(withCASPT2) inputFile.println("&caspt2");
				inputFile.close();
				return true;
			} else if(stateNames[state].trim().equals("D+")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
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
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				if(withCASPT2) inputFile.println("&caspt2");
				inputFile.close();
				return true;
			} else if(stateNames[state].trim().equals("T+")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				inputFile.println("nactel");
				inputFile.println(" "+(numCASe-1));
				inputFile.println("spin");
				inputFile.println(" 3");
				inputFile.println("inactive");
				inputFile.println(" "+Inact);
				inputFile.println("ras2");
				inputFile.println(" "+numCASo);
				inputFile.println("prwf");
				inputFile.println("  0");
				inputFile.println("prsd");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				if(withCASPT2) inputFile.println("&caspt2");
				inputFile.close();
				return true;
			} else if(stateNames[state].trim().equals("S-")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				inputFile.println("nactel");
				inputFile.println(" "+(numCASe+1));
				inputFile.println("spin");
				inputFile.println(" 1");
				inputFile.println("inactive");
				inputFile.println(" "+Inact);
				inputFile.println("ras2");
				inputFile.println(" "+numCASo);
				inputFile.println("prwf");
				inputFile.println("  0");
				inputFile.println("prsd");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				if(withCASPT2) inputFile.println("&caspt2");
				inputFile.close();
				return true;
			} else if(stateNames[state].trim().equals("D-")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
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
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				if(withCASPT2) inputFile.println("&caspt2");
				inputFile.close();
				return true;
			} else if(stateNames[state].trim().equals("T-")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				inputFile.println("nactel");
				inputFile.println(" "+(numCASe+1));
				inputFile.println("spin");
				inputFile.println(" 3");
				inputFile.println("inactive");
				inputFile.println(" "+Inact);
				inputFile.println("ras2");
				inputFile.println(" "+numCASo);
				inputFile.println("prwf");
				inputFile.println("  0");
				inputFile.println("prsd");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				if(withCASPT2) inputFile.println("&caspt2");
				inputFile.close();
				return true;
			}
		} catch(IOException e) {
			return false;
		}
		return false;
	}
	
	public Boolean Molcas_SCF_Converged(Integer frag, Integer numCASe) {
		String fileName = projectRoot+fragmentNames[frag]+"_SCF.output";
		String card;
		Integer numOcc;
		Boolean converged = false;
		Double energy = 0.0;
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
	
			numOcc=0;
			Integer numSec=0;
			Integer numBas=0;
			Integer numOrb=0;
			Integer index=0;
			while((card=br.readLine()) != null) {
				if(card.contains("      Occupied orbitals")) {
					numOcc=Integer.valueOf(card.substring(33,38).trim());
				}
				if(card.contains("      Secondary orbitals")) {
					numSec=Integer.valueOf(card.substring(33,38).trim());
				}
				if(card.contains("::    Total SCF energy")) {
					converged=true;
					energy=Double.valueOf(card.substring(51,68).trim()).doubleValue();
				}
				if(card.contains("++    Molecular orbitals:") && converged) {
					numBas=numOcc+numSec;
					numOrb=2*numOcc;
					Double[] occ = new Double[numOrb];
					Integer[] typ = new Integer[numOrb];
					Double[] coef = new Double[10];
					Double[] cmax = new Double[10];
					Integer[] imax = new Integer[10];
					Integer ilast;
					Integer count=0;
					card=br.readLine();
					card=br.readLine();
					card=br.readLine();
					card=br.readLine();
					card=br.readLine();
					card=br.readLine();
					card=br.readLine();
					card=br.readLine();
					card=br.readLine();
					while(count<numOrb) {
						ilast=Integer.min(10,numOrb-count);
						for(int i=0; i<ilast; i++) {
							occ[count+i]=Double.valueOf(card.substring(24+i*10,34+i*10).trim()).doubleValue();
							typ[count+i]=0;
							cmax[i]=0.0;
						}
						card=br.readLine();
						for(int j=0; j<numBas; j++) {
							card=br.readLine();
							for(int i=0; i<ilast; i++) {
								coef[i]=Double.valueOf(card.substring(24+i*10,34+i*10).trim()).doubleValue();
								if(Math.abs(coef[i])>cmax[i]) {
									imax[i]=j;
									if(card.contains("s")) typ[count+i]=1;
									if(card.contains("px")) typ[count+i]=2;
									if(card.contains("py")) typ[count+i]=3;
									if(card.contains("pz")) typ[count+i]=4;
									if(card.contains("d2-")) typ[count+i]=5;
									if(card.contains("d1-")) typ[count+i]=6;
									if(card.contains("d0")) typ[count+i]=7;
									if(card.contains("d1+")) typ[count+i]=8;
									if(card.contains("d2+")) typ[count+i]=9;
								}
								cmax[i]=Double.max(cmax[i],Math.abs(coef[i]));
							}	
						}
						card=br.readLine();
						card=br.readLine();
						card=br.readLine();
						card=br.readLine();
						count=count+10;
					}
					Integer n=numCASe/2, numpz=0;
					for(int i=numOcc-n; i<numOcc; i++) {
						index=i;
						if(typ[index]==4) numpz++;
					}
					numAlt=0;
					if(numpz==n) {
						for(int i=numOcc; i<numOcc+n; i++) {
							index=i;
							if(typ[index]!=4) {
								alter[numAlt][0]=i;
								numAlt++;
							}
						}
						Integer ii=numOcc+n;
						for(int ialt=0; ialt<numAlt; ialt++) {
							while(typ[ii]!=4) ii++;
							alter[ialt][1]=ii; ii++;
						}
					}
				}
			}
			br.close();
			return converged;
		} catch(IOException ef) {
			return false;
		}
	}

	public Double Molcas_SCF(Integer frag, Integer numCASe) {
		String fileName = projectRoot+fragmentNames[frag]+"_SCF.output";
		String card;
		Integer numOcc;
		Boolean converged = false;
		Double energy = 0.0;
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
	
			numOcc=0;
			Integer numSec=0;
			Integer numBas=0;
			Integer numOrb=0;
			Integer index=0;
			while((card=br.readLine()) != null) {
				if(card.contains("      Occupied orbitals")) {
					numOcc=Integer.valueOf(card.substring(33,38).trim());
				}
				if(card.contains("      Secondary orbitals")) {
					numSec=Integer.valueOf(card.substring(33,38).trim());
				}
				if(card.contains("::    Total SCF energy")) {
					converged=true;
					energy=Double.valueOf(card.substring(51,68).trim()).doubleValue();
				}
				if(card.contains("++    Molecular orbitals:") && converged) {
					numBas=numOcc+numSec;
					numOrb=2*numOcc;
					Double[] occ = new Double[numOrb];
					Integer[] typ = new Integer[numOrb];
					Double[] coef = new Double[10];
					Double[] cmax = new Double[10];
					Integer[] imax = new Integer[10];
					Integer ilast;
					Integer count=0;
					card=br.readLine();
					card=br.readLine();
					card=br.readLine();
					card=br.readLine();
					card=br.readLine();
					card=br.readLine();
					card=br.readLine();
					card=br.readLine();
					card=br.readLine();
					while(count<numOrb) {
						ilast=Integer.min(10,numOrb-count);
						for(int i=0; i<ilast; i++) {
							occ[count+i]=Double.valueOf(card.substring(24+i*10,34+i*10).trim()).doubleValue();
							typ[count+i]=0;
							cmax[i]=0.0;
						}
						card=br.readLine();
						for(int j=0; j<numBas; j++) {
							card=br.readLine();
							for(int i=0; i<ilast; i++) {
								coef[i]=Double.valueOf(card.substring(24+i*10,34+i*10).trim()).doubleValue();
								if(Math.abs(coef[i])>cmax[i]) {
									imax[i]=j;
									if(card.contains("s")) typ[count+i]=1;
									if(card.contains("px")) typ[count+i]=2;
									if(card.contains("py")) typ[count+i]=3;
									if(card.contains("pz")) typ[count+i]=4;
									if(card.contains("d2-")) typ[count+i]=5;
									if(card.contains("d1-")) typ[count+i]=6;
									if(card.contains("d0")) typ[count+i]=7;
									if(card.contains("d1+")) typ[count+i]=8;
									if(card.contains("d2+")) typ[count+i]=9;
								}
								cmax[i]=Double.max(cmax[i],Math.abs(coef[i]));
							}	
						}
						card=br.readLine();
						card=br.readLine();
						card=br.readLine();
						card=br.readLine();
						count=count+10;
					}
					Integer n=numCASe/2, numpz=0;
					for(int i=numOcc-n; i<numOcc; i++) {
						index=i;
						if(typ[index]==4) numpz++;
					}
					numAlt=0;
					if(numpz==n) {
						for(int i=numOcc; i<numOcc+n; i++) {
							index=i;
							if(typ[index]!=4) {
								alter[numAlt][0]=i;
								numAlt++;
							}
						}
						Integer ii=numOcc+n;
						for(int ialt=0; ialt<numAlt; ialt++) {
							while(typ[ii]!=4) ii++;
							alter[ialt][1]=ii; ii++;
						}
					}
				}
			}
			br.close();
			return energy;
		} catch(IOException ef) {
			return energy;
		}
	}
	
	public Integer Molcas_CASSCF_Converged(Integer frag, Integer state) {
		String fileName = projectRoot+fragmentNames[frag]+"_"+stateNames[state]+".output";
		String card;
		Boolean convergedCASSCF = false;
		Integer numConverged = 0;
		Double energyCASSCF = 0.0;
		Double energyCASPT2 = 0.0;
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			while((card=br.readLine()) != null) {
				if(card.startsWith("      RASSCF energy for state")) {
					numConverged=1;
					energyCASSCF=Double.valueOf(card.substring(48,64).trim()).doubleValue();
				}
				if(convergedCASSCF && card.startsWith("      Total energy:    ")) {
					numConverged=2;
					energyCASPT2=Double.valueOf(card.substring(28,46).trim()).doubleValue();
				}
			}
			br.close();
			return numConverged;
		} catch(IOException ef) {
			return 0;
		}
	}

	public Double Molcas_CASSCF(Integer frag, Integer state) {
		String fileName = projectRoot+fragmentNames[frag]+"_"+stateNames[state]+".output";
		String card;
		Boolean convergedCASSCF = false;
		Integer numConverged = 0;
		Double energyCASSCF = 0.0;
		Double energyCASPT2 = 0.0;
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			while((card=br.readLine()) != null) {
				if(card.startsWith("      RASSCF energy for state")) {
					numConverged=1;
					energyCASSCF=Double.valueOf(card.substring(48,64).trim()).doubleValue();
				}
				if(numConverged==1 && card.startsWith("      Total energy:    ")) {
					numConverged=2;
					energyCASPT2=Double.valueOf(card.substring(28,46).trim()).doubleValue();
				}
			}
			br.close();
			return energyCASSCF;
		} catch(IOException ef) {
			return 0.0;
		}
	}
	
	public Double Molcas_CASPT2(Integer frag, Integer state) {
		String fileName = projectRoot+fragmentNames[frag]+"_"+stateNames[state]+".output";
		String card;
		Boolean convergedCASSCF = false;
		Integer numConverged = 0;
		Double energyCASSCF = 0.0;
		Double energyCASPT2 = 0.0;
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			while((card=br.readLine()) != null) {
				if(card.startsWith("      RASSCF energy for state")) {
					numConverged=1;
					energyCASSCF=Double.valueOf(card.substring(48,64).trim()).doubleValue();
				}
				if(numConverged==1 && card.startsWith("      Total energy:    ")) {
					numConverged=2;
					energyCASPT2=Double.valueOf(card.substring(28,46).trim()).doubleValue();
				}
			}
			br.close();
			return energyCASPT2;
		} catch(IOException ef) {
			return 0.0;
		}
	}

	public Boolean write_MEBF_XYZ(String p, String pn, Integer n, String[] frags) {
		String fileName = p+".xyz";
		Integer numberAtoms = 0;
		for(int j=0; j<n; j++) {
			initialize2(pn, frags[j], 2);
			numberAtoms=numberAtoms+numAtoms;
		}
		try {
			PrintfWriter xyzFile = new PrintfWriter(new FileWriter(fileName));
			xyzFile.printf("%6d",numberAtoms);
			xyzFile.println();
			xyzFile.println("Coordinates in Angstrom");
			for(int j=0; j<n; j++) {
				initialize2(pn, frags[j], 2);
				for(int i=0; i<numAtoms; i++) {
					xyzFile.print(atomLabel[i]); 
					xyzFile.printf("%12.8f",coordinates[i][0]);
					xyzFile.printf("%12.8f",coordinates[i][1]);
					xyzFile.printf("%12.8f",coordinates[i][2]);
					xyzFile.println();
				}
			}
			xyzFile.close();
			return true;
		} catch(IOException ei) {
			return false;
		}
	}
	
	public void write_Molcas_MEBF_One(String p, String pn, Integer n, String[] frags) {
		String fileName = p+"_MEBFONE.input";
		String previous;
		Integer atomNumber=0;
		Integer maxElement = 45;
		Integer[] counts = new Integer[maxElement];
		for(int i=0; i<maxElement; i++) counts[i]=0;
		try {
			PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
			inputFile.println("&seward");
			inputFile.println("high cholesky");
			for(int j=0; j<n; j++) {
				initialize2(pn, frags[j], 2);
			    previous=" ";
			    for(int i=0; i<numAtoms; i++) {
			    	if(!atomLabel[i].trim().equals(previous)) {
			    		previous=atomLabel[i].trim();
			    		if(i>0) inputFile.println("end basis set");
			    		inputFile.println("basis set");
			    		if(previous.equals("H")) {
			    			inputFile.println("h.ano-s...3s3p.");
			    			atomNumber=1;
			    		}
			    		if(previous.equals("C")) {
			    			inputFile.println("c.ano-s...4s3p2d.");
			    			atomNumber=6;
			    		}
			    		if(previous.equals("N")) {
			    			inputFile.println("n.ano-s...4s3p2d.");
			    			atomNumber=7;
			    		}
			    		if(previous.equals("O")) {
			    			inputFile.println("o.ano-s...4s3p2d.");
			    			atomNumber=8;
			    		}
			    		if(previous.equals("F")) {
			    			inputFile.println("f.ano-s...4s3p2d.");
			    			atomNumber=9;
			    		}
			    		if(previous.equals("Cl")) {
			    			inputFile.println("cl.ano-s...4s3p2d.");
			    			atomNumber=17;
			    		}
			    		if(previous.equals("Br")) {
			    			inputFile.println("br.ano-s...4s3p2d.");
			    			atomNumber=35;
			    		}
			    	}
			    	counts[atomNumber]++;
			    	inputFile.print(atomLabel[i].trim()+counts[atomNumber]+" "); 
					inputFile.printf("%12.8f",coordinates[i][0]);
					inputFile.printf("%12.8f",coordinates[i][1]);
					inputFile.printf("%12.8f",coordinates[i][2]);
					inputFile.println(" Angstrom");
					previous=atomLabel[i].trim();
			    }
			    inputFile.println("end basis set");
			}
		    inputFile.println("oneonly");
		    inputFile.println(">>> COPY "+p+".RunFile $CurrDir/"+p+".RUNFILE");
		    inputFile.close();
//		    return true;
		} catch(IOException e) {
//			return false;
		}
	}

	public void write_Molcas_MEBF_Two(String p, String pn, Integer n, String[] frags) {
		String fileName = p+"_MEBFTWO.input";
		String previous;
		Integer atomNumber=0;
		Integer maxElement = 45;
		Integer[] counts = new Integer[maxElement];
		for(int i=0; i<maxElement; i++) counts[i]=0;
		try {
			PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
			inputFile.println("&seward");
			inputFile.println("high cholesky");
			for(int j=0; j<n; j++) {
				initialize2(pn, frags[j], 2);
			    previous=" ";
			    for(int i=0; i<numAtoms; i++) {
			    	if(!atomLabel[i].trim().equals(previous)) {
			    		previous=atomLabel[i].trim();
			    		if(i>0) inputFile.println("end basis set");
			    		inputFile.println("basis set");
			    		if(previous.equals("H")) {
			    			inputFile.println("h.ano-s...3s3p.");
			    			atomNumber=1;
			    		}
			    		if(previous.equals("C")) {
			    			inputFile.println("c.ano-s...4s3p2d.");
			    			atomNumber=6;
			    		}
			    		if(previous.equals("N")) {
			    			inputFile.println("n.ano-s...4s3p2d.");
			    			atomNumber=7;
			    		}
			    		if(previous.equals("O")) {
			    			inputFile.println("o.ano-s...4s3p2d.");
			    			atomNumber=8;
			    		}
			    		if(previous.equals("F")) {
			    			inputFile.println("f.ano-s...4s3p2d.");
			    			atomNumber=9;
			    		}
			    		if(previous.equals("Cl")) {
			    			inputFile.println("cl.ano-s...4s3p2d.");
			    			atomNumber=17;
			    		}
			    		if(previous.equals("Br")) {
			    			inputFile.println("br.ano-s...4s3p2d.");
			    			atomNumber=35;
			    		}
			    	}
			    	counts[atomNumber]++;
			    	inputFile.print(atomLabel[i].trim()+counts[atomNumber]+" "); 
					inputFile.printf("%12.8f",coordinates[i][0]);
					inputFile.printf("%12.8f",coordinates[i][1]);
					inputFile.printf("%12.8f",coordinates[i][2]);
					inputFile.println(" Angstrom");
					previous=atomLabel[i].trim();
			    }
			    inputFile.println("end basis set");   
			}
		    inputFile.println(">>> COPY $CurrDir/COMMONORB INPORB");
		    inputFile.println(">>> COPY $CurrDir/COMMONORB $CurrDir/"+p+".COMMONORB");
		    inputFile.println();
		    inputFile.println("&motra");
		    inputFile.println("LumOrb");
		    inputFile.println("frozen");
		    inputFile.println(" 0");
		    inputFile.println("deleted");
		    inputFile.println(" $DELETED");
		    inputFile.println("ctonly");
		    inputFile.println("kpq");
		    inputFile.println();
		    inputFile.println(">>> COPY "+p+".RunFile $CurrDir/"+p+".RUNFILE");
		    inputFile.println(">>> COPY "+p+".OneInt  $CurrDir/"+p+".ONEINT");
		    inputFile.println(">>> COPY "+p+".TraOne  $CurrDir/"+p+".TRAONE");
		    inputFile.println(">>> COPY "+p+".ChVec1  $CurrDir/"+p+".CHVEC1");
		    inputFile.println(">>> COPY "+p+".ChRed   $CurrDir/"+p+".CHRED");
		    inputFile.println(">>> COPY "+p+".ChRst   $CurrDir/"+p+".CHORST");
		    inputFile.println(">>> COPY "+p+".ChMap   $CurrDir/"+p+".CHOMAP");
		    inputFile.println(">>> COPY _CHMOT1        $CurrDir/"+p+"_CHMOT1");
		    inputFile.println(">>> eval NPROCS = $MOLCAS_NPROCS - 1");
		    inputFile.println(">>> foreach L in (1 .. $NPROCS )");
		    inputFile.println(">>> shell cat tmp_$L/_CHMOT1 >> $CurrDir/"+p+"_CHMOT1");
		    inputFile.println(">>> enddo");
		    inputFile.close();
//		    return true;
		} catch(IOException e) {
//			return false;
		}
	}

	public void write_Molcas_MEBF_CB(String p, Integer n, Integer[] states) {
		String fileName = p+"_MEBFCB.input";
		try {
			PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
			inputFile.println("Project");
			inputFile.println(p);
			inputFile.println("Fragments");
			inputFile.printf("%3d",n); inputFile.println();
			for(int i=0; i<n; i++) inputFile.printf("%3d",states[i]); inputFile.println();
			inputFile.println("Threshold");
			inputFile.println(" 1.0e-5");
			inputFile.println("Labels");
			for(int i=0; i<n; i++) {
				for(int j=0; j<states[i]; j++) inputFile.print(" "+stateNames[j].trim());
			}
			inputFile.println();
			inputFile.println("Energies");
		    inputFile.close();
			} catch(IOException e) {
			}
		
	}

	public void rotate_AND_translate(Double[] rt) {
		Double tx=rt[0];
		Double ty=rt[1];
		Double tz=rt[2];
		Double rx=rt[3];
		Double ry=rt[4];
		Double rz=rt[5];
		v[0]=0.0;
		v[1]=0.0;
		v[2]=0.0;
		w[0]=0.0;
		w[1]=0.0;
		w[2]=1.0;
		g=rz;
		for(int i=0; i<numAtoms; i++) {
			x[0]=coordinates[i][0];
			x[1]=coordinates[i][1];
			x[2]=coordinates[i][2];
			rotate();
			coordinates[i][0]=y[0];
			coordinates[i][1]=y[1];
			coordinates[i][2]=y[2];
		}
		w[0]=0.0;
		w[1]=1.0;
		w[2]=0.0;
		g=ry;
		for(int i=0; i<numAtoms; i++) {
			x[0]=coordinates[i][0];
			x[1]=coordinates[i][1];
			x[2]=coordinates[i][2];
			rotate();
			coordinates[i][0]=y[0];
			coordinates[i][1]=y[1];
			coordinates[i][2]=y[2];
		}
		w[0]=1.0;
		w[1]=0.0;
		w[2]=0.0;
		g=rx;
		for(int i=0; i<numAtoms; i++) {
			x[0]=coordinates[i][0];
			x[1]=coordinates[i][1];
			x[2]=coordinates[i][2];
			rotate();
			coordinates[i][0]=y[0];
			coordinates[i][1]=y[1];
			coordinates[i][2]=y[2];
		}
		for(int i=0; i<numAtoms; i++) {
			coordinates[i][0]=coordinates[i][0]+tx;
			coordinates[i][1]=coordinates[i][1]+ty;
			coordinates[i][2]=coordinates[i][2]+tz;
		}
	}
	
	public void rotate() {
		Double small = 1.0e-24;
		Double a, b, pi, r;
		Double sa, ca, sb, cb, sg, cg;
		Double[] xx = new Double[3];
		Double[][] t = new Double[3][3];
		if(Math.abs(w[1]-v[1])<small) {
			if(Math.abs(w[0]-v[0])<small) {
				a = 0.0;
			} else {
				if(w[0]-v[0]>0.0) {
					a=2.0*Math.atan(1.0);
				} else {
					a=-2.0*Math.atan(1.0);
				}
			}
		} else {
			a = Math.atan(Math.abs(w[0]-v[0])/Math.abs(w[1]-v[1]));
			pi = 4.0*Math.atan(1.0);
			if(w[0]-v[0]>0.0 && w[1]-v[1]<0.0) a=pi-a;
			if(w[0]-v[0]<0.0 && w[1]-v[1]>0.0) a=-a;
			if(w[0]-v[0]<0.0 && w[1]-v[1]<0.0) a=pi+a;
		}
		r=(w[0]-v[0])*(w[0]-v[0])+(w[1]-v[1])*(w[1]-v[1])+(w[2]-v[2])*(w[2]-v[2]);
		for(int i=0; i<3; i++) xx[i]=x[i]-v[i];
		if(r<small) {
			for(int i=0; i<3; i++) y[i]=x[i];
		} else {
			b=Math.acos((w[2]-v[2])/Math.sqrt(r));
			sa=Math.sin(a);
			ca=Math.cos(a);
			sb=Math.sin(b);
			cb=Math.cos(b);
			sg=Math.sin(g);
			cg=Math.cos(g);
	        t[0][0]=ca*ca*cg-sa*ca*cb*sg+sa*ca*cb*sg+sa*sa*cb*cb*cg+sa*sa*sb*sb;
	        t[0][1]=(-sa)*ca*cg-ca*ca*cb*sg-sa*sa*cb*sg+sa*ca*cb*cb*cg+sa*ca*sb*sb;
	        t[0][2]=ca*sb*sg-sa*sb*cb*cg+sa*sb*cb;
	        t[1][0]=(-sa)*ca*cg+sa*sa*cb*sg+ca*ca*cb*sg+sa*ca*cb*cb*cg+sa*ca*sb*sb;
	        t[1][1]=sa*sa*cg+sa*ca*cb*sg-sa*ca*cb*sg+ca*ca*cb*cb*cg+ca*ca*sb*sb;
	        t[1][2]=(-sa)*sb*sg-ca*sb*cb*cg+ca*sb*cb;
	        t[2][0]=(-ca)*sb*sg-sa*sb*cb*cg+sa*sb*cb;
	        t[2][1]=sa*sb*sg-ca*sb*cb*cg+ca*sb*cb;
	        t[2][2]=sb*sb*cg+cb*cb;
	        for(int i=0; i<3; i++) y[i]=xx[0]*t[i][0]+xx[1]*t[i][1]+xx[2]*t[i][2]+v[i];
		}
	}

}