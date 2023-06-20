package gronor;

import java.util.*;
import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;

import cformat.PrintfWriter;

public class gronor_Fragment {

	private static final long serialVersionUID = 4L;

	String projectName;
	
	Integer maxAtoms = 100;
	
	String projectRoot;
	
	String[] fragmentNames = new String[] {"A", "B", "C", "D", "E", "F", "G", "H"};
    String[] stateNames = new String[] {"S0","S1","S2","D0","D1","T1","T2","S+","D+","T+","S-","D-","T-","q1","Q1","SQ1"};
    
    
    String[]  elements = new String[] {"H", "C", "N", "O", "F", "Cl", "Cu", "Se", "Br"};
    Integer[] atomNumber = new Integer[] {1, 6, 7, 8, 9, 17, 29, 34, 35};
    String[] basisSet = new String[] {" ", " ", " ", " ", " ", " ", " ", " ", " "};
	
	Integer fragmentID = 0, numAtoms=0, numOcc=0;
	String fragmentName = " ";
	String fragmentRoot = " ";
	String fragmentAXYZ, fragmentXYZ, nwchemOutput;

	String convergedName = " ";
	
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
	Boolean altDone=false;

	gronor_Fragment(){
	}

	public void initialize(String nameP, String nameF,String nameA,String nameB, Integer numE, Double[] rt) {
		// Initialize PA from source FB
		Boolean writeFile = true; 
		projectRoot=nameP;
		fragmentRoot=nameF;
		int last=nameF.indexOf("_");
		if(last>0) fragmentRoot=nameF.substring(0,last);

		if(fragmentRoot.length()==0) return;
		
		fragmentName=nameF.trim();
		numOcc = numE / 2;
		for(int i=0; i<6; i++) RandT[i]=rt[i];
			if(!nameF.trim().equals("")) {
				if(read_XYZ()) {
					rotate_AND_translate(RandT);
					fragmentName=nameP.trim()+nameA.trim();
					if(writeFile) {
						if(!nameA.trim().equals("")) write_XYZ();
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
	}

	public void initialize2(String nameF, String nameA, Integer numE) {
		fragmentRoot=nameF;
		fragmentName=projectRoot.trim()+nameA.trim();
		numOcc = numE / 2;
		if(read_XYZ()) {
			rotate_AND_translate(RandT);
		} else {
		}
	}

	public void initialize3(String nameF, String nameA, Integer numE) {
		fragmentRoot=nameF;
		fragmentName=projectRoot.trim()+nameA.trim();
		numOcc = numE / 2;
		if(read_XYZ()) {
		} else {
			System.out.println("XYZ file "+fragmentName+".xyz not found"); 
			System.exit(0);
		}
	}

	public void write_Run_Script_Fragments(String fn, String fp, String fs, Integer frag, Integer stateset, Integer[] lenStateList, Integer[][] ndxStateList, Integer ranks, Integer memory, String acc, String jobName, String limit) {
		String rootName = projectRoot.trim()+fragmentNames[frag].trim();
		String fileName = fp.trim()+fs.trim()+"_SCF.run";
		String slurmName = fp.trim()+fs.trim()+"_SCF.slurm";
		String lsfName = fp.trim()+fs.trim()+"_SCF.lsf";
		String fullName;
		Integer numStates = lenStateList[stateset];
		try {
			PrintfWriter runFile = new PrintfWriter(new FileWriter(fileName));
			PrintfWriter slurmFile = new PrintfWriter(new FileWriter(slurmName));
			PrintfWriter lsfFile = new PrintfWriter(new FileWriter(lsfName));
			runFile.println("#!/usr/bin/tcsh");
			runFile.println("setenv MOLCAS_NPROCS "+ranks);
			runFile.println("setenv MOLCAS_MEM "+memory);

			slurmFile.println("#!/usr/bin/tcsh");
			slurmFile.println("#SBATCH --account="+acc.trim());
			slurmFile.println("#SBATCH --nodes=1");
			slurmFile.println("#SBATCH --tasks-per-node="+ranks);
			slurmFile.println("#SBATCH --job-name="+jobName.trim());
			slurmFile.println("#SBATCH --time="+limit.trim());
			slurmFile.println();
			slurmFile.println("module purge");
			slurmFile.println("module load intel openmpi");
			slurmFile.println();
			slurmFile.println("#source required files");
			slurmFile.println();
			slurmFile.println("#cd to job directory");
			slurmFile.println();

			lsfFile.println("#!/usr/bin/tcsh");
			lsfFile.println("#BSUB -P "+acc.trim());
			lsfFile.println("#BSUB -J "+jobName.trim());
			lsfFile.println("#BSUB -W "+limit.trim());
			lsfFile.println("#BSUB -q batch ");
			lsfFile.println("#BSUB -o job.o%J");
			lsfFile.println("#BSUB -e job.e%J");
			lsfFile.println("#BSUB -alloc_flags \"gpumps\"");
			lsfFile.println("#BSUB -nnodes "+ranks);
			lsfFile.println();
			lsfFile.println("#module purge");
			lsfFile.println("#module load intel openmpi");
			lsfFile.println();
			lsfFile.println("#source required files");
			lsfFile.println();
			lsfFile.println("#cd to job directory");
			lsfFile.println();
			
			fullName = fp.trim()+fs.trim()+"_INT";
			runFile.println("cp "+fullName.trim()+".input "+fp.trim()+fs.trim()+".input; "+"pymolcas "+fp.trim()+fs.trim()+".input > "+fullName.trim()+".output");
			slurmFile.println("cp "+fullName.trim()+".input "+fp.trim()+fs.trim()+".input; "+"pymolcas "+fp.trim()+fs.trim()+".input > "+fullName.trim()+".output");
			lsfFile.println("cp "+fullName.trim()+".input "+fp.trim()+fs.trim()+".input; "+"pymolcas "+fp.trim()+fs.trim()+".input > "+fullName.trim()+".output");
			fullName = fp.trim()+fs.trim()+"_SCF";
			runFile.println("cp "+fullName.trim()+".input "+fp.trim()+fs.trim()+".input; "+"pymolcas "+fp.trim()+fs.trim()+".input > "+fullName.trim()+".output");
			slurmFile.println("cp "+fullName.trim()+".input "+fp.trim()+fs.trim()+".input; "+"pymolcas "+fp.trim()+fs.trim()+".input > "+fullName.trim()+".output");
			lsfFile.println("cp "+fullName.trim()+".input "+fp.trim()+fs.trim()+".input; "+"pymolcas "+fp.trim()+fs.trim()+".input > "+fullName.trim()+".output");
			runFile.println("alter "+fullName.trim()+".output "+fp.trim()+fs.trim()+".alter");
			slurmFile.println("alter "+fullName.trim()+".output "+fp.trim()+fs.trim()+".alter");
			lsfFile.println("alter "+fullName.trim()+".output "+fp.trim()+fs.trim()+".alter");
			
			runFile.println("rm "+fp.trim()+fs.trim()+".input");
			slurmFile.println("rm "+fp.trim()+fs.trim()+".input");
			lsfFile.println("rm "+fp.trim()+fs.trim()+".input");
			
			runFile.close();
			slurmFile.close();
			lsfFile.close();
		} catch(IOException ei) {
		}
		fileName = fp.trim()+fs.trim()+".run";
		slurmName = fp.trim()+fs.trim()+".slurm";
		lsfName = fp.trim()+fs.trim()+".lsf";
		try {
			PrintfWriter runFile = new PrintfWriter(new FileWriter(fileName));
			PrintfWriter slurmFile = new PrintfWriter(new FileWriter(slurmName));
			PrintfWriter lsfFile = new PrintfWriter(new FileWriter(lsfName));
			runFile.println("#!/usr/bin/tcsh");
			runFile.println("setenv MOLCAS_NPROCS "+ranks);
			runFile.println("setenv MOLCAS_MEM "+memory);
			
			slurmFile.println("#!/usr/bin/tcsh");
			slurmFile.println("#SBATCH --account="+acc.trim());
			slurmFile.println("#SBATCH --nodes=1");
			slurmFile.println("#SBATCH --tasks-per-node="+ranks);
			slurmFile.println("#SBATCH --job-name="+jobName.trim());
			slurmFile.println("#SBATCH --time="+limit.trim());
			slurmFile.println();
			slurmFile.println("module purge");
			slurmFile.println("module load intel openmpi");
			slurmFile.println();
			slurmFile.println("#source required files");
			slurmFile.println();
			slurmFile.println("#cd to job directory");
			slurmFile.println();

			lsfFile.println("#!/usr/bin/tcsh");
			lsfFile.println("#BSUB -P "+acc.trim());
			lsfFile.println("#BSUB -J "+jobName.trim());
			lsfFile.println("#BSUB -W "+limit.trim());
			lsfFile.println("#BSUB -q batch ");
			lsfFile.println("#BSUB -o job.o%J");
			lsfFile.println("#BSUB -e job.e%J");
			lsfFile.println("#BSUB -alloc_flags \"gpumps\"");
			lsfFile.println("#BSUB -nnodes "+ranks);
			lsfFile.println();
			lsfFile.println("#module purge");
			lsfFile.println("#module load intel openmpi");
			lsfFile.println();
			lsfFile.println("#source required files");
			lsfFile.println();
			lsfFile.println("#cd to job directory");
			lsfFile.println();
			
			runFile.println("# OpenMolcas runs for the spin states are run in the same workspace "+fp.trim()+fs.trim());
			slurmFile.println("# OpenMolcas runs for the spin states are run in the same workspace "+fp.trim()+fs.trim());
			lsfFile.println("# OpenMolcas runs for the spin states are run in the same workspace "+fp.trim()+fs.trim());
			
			runFile.println("# Evaluate the integrals");
			slurmFile.println("# Evaluate the integrals");
			lsfFile.println("# Evaluate the integrals");
			fullName = fp.trim()+fs.trim()+"_INT";
			runFile.println("cp "+fullName.trim()+".input "+fp.trim()+fs.trim()+".input; "+"pymolcas "+fp.trim()+fs.trim()+".input > "+fullName.trim()+".output");
			slurmFile.println("cp "+fullName.trim()+".input "+fp.trim()+fs.trim()+".input; "+"pymolcas "+fp.trim()+fs.trim()+".input > "+fullName.trim()+".output");
			lsfFile.println("cp "+fullName.trim()+".input "+fp.trim()+fs.trim()+".input; "+"pymolcas "+fp.trim()+fs.trim()+".input > "+fullName.trim()+".output");
			runFile.println("# SCF to generate starting orbitals");
			slurmFile.println("# SCF to generate starting orbitals");
			lsfFile.println("# SCF to generate starting orbitals");
			fullName = fp.trim()+fs.trim()+"_SCF";
			runFile.println("cp "+fullName.trim()+".input "+fp.trim()+fs.trim()+".input; "+"pymolcas "+fp.trim()+fs.trim()+".input > "+fullName.trim()+".output");
			slurmFile.println("cp "+fullName.trim()+".input "+fp.trim()+fs.trim()+".input; "+"pymolcas "+fp.trim()+fs.trim()+".input > "+fullName.trim()+".output");
			lsfFile.println("cp "+fullName.trim()+".input "+fp.trim()+fs.trim()+".input; "+"pymolcas "+fp.trim()+fs.trim()+".input > "+fullName.trim()+".output");
			runFile.println("# RASSCF and CASPT2 for each spin state");
			slurmFile.println("# RASSCF and CASPT2 for each spin state");
			lsfFile.println("# RASSCF and CASPT2 for each spin state");
			for(int i=0; i<numStates; i++) {
				Integer index = ndxStateList[stateset][i];
				fullName = fp.trim()+fs.trim()+"_"+stateNames[index].trim();
				runFile.println("cp "+fullName.trim()+".input "+fp.trim()+fs.trim()+".input; "+"pymolcas "+fp.trim()+fs.trim()+".input > "+fullName.trim()+".output");
				slurmFile.println("cp "+fullName.trim()+".input "+fp.trim()+fs.trim()+".input; "+"pymolcas "+fp.trim()+fs.trim()+".input > "+fullName.trim()+".output");
				lsfFile.println("cp "+fullName.trim()+".input "+fp.trim()+fs.trim()+".input; "+"pymolcas "+fp.trim()+fs.trim()+".input > "+fullName.trim()+".output");
			}
			runFile.println("# Remove the generic input file to avoid clutter");
			slurmFile.println("# Remove the generic input file to avoid clutter");
			lsfFile.println("# Remove the generic input file to avoid clutter");
			runFile.println("rm "+fp.trim()+fs.trim()+".input");
			slurmFile.println("rm "+fp.trim()+fs.trim()+".input");
			lsfFile.println("rm "+fp.trim()+fs.trim()+".input");
			runFile.println("# Parse the output files for RASSCF and CASPT2 energies and add them to the det file headers");
			slurmFile.println("# Parse the output files for RASSCF and CASPT2 energies and add them to the det file headers");
			lsfFile.println("# Parse the output files for RASSCF and CASPT2 energies and add them to the det file headers");
			for(int i=0; i<numStates; i++) {
				Integer index = ndxStateList[stateset][i];
				fullName = fp.trim()+fs.trim()+"_"+stateNames[index].trim();
				runFile.println("det_header "+fullName.trim()+".output "+fullName.trim()+".det");
				slurmFile.println("det_header "+fullName.trim()+".output "+fullName.trim()+".det");
				lsfFile.println("det_header "+fullName.trim()+".output "+fullName.trim()+".det");
			}
			
			runFile.close();
			slurmFile.close();
			lsfFile.close();
		} catch(IOException ei) {
		}
	}

	public void write_Rotate_Script_Fragments(String fn, String fp, String fs, Integer frag, Integer stateset, Integer[] lenStateList, Integer[][] ndxStateList, Object name1, Object name2, Object nameF, Integer ranks, Integer memory, String acc, String jobName, String limit) {
		String rootName = projectRoot.trim()+fragmentNames[frag].trim();
		Integer numStates = lenStateList[stateset];
		String fileName = fp.trim()+fs.trim()+".run";
		String slurmName = fp.trim()+fs.trim()+".slurm";
		String lsfName = fp.trim()+fs.trim()+".lsf";
		String fullName;
		try {
			PrintfWriter runFile = new PrintfWriter(new FileWriter(fileName));
			PrintfWriter slurmFile = new PrintfWriter(new FileWriter(slurmName));
			PrintfWriter lsfFile = new PrintfWriter(new FileWriter(lsfName));
			runFile.println("#!/usr/bin/tcsh");
			runFile.println("setenv MOLCAS_NPROCS "+ranks);
			runFile.println("setenv MOLCAS_MEM "+memory);
			
			slurmFile.println("#!/usr/bin/tcsh");
			slurmFile.println("#SBATCH --account="+acc.trim());
			slurmFile.println("#SBATCH --nodes=1");
			slurmFile.println("#SBATCH --tasks-per-node="+ranks);
			slurmFile.println("#SBATCH --job-name="+jobName.trim());
			slurmFile.println("#SBATCH --time="+limit.trim());
			slurmFile.println();
			slurmFile.println("module purge");
			slurmFile.println("module load intel openmpi");
			slurmFile.println();
			slurmFile.println("#source required files");
			slurmFile.println();
			slurmFile.println("#cd to job directory");
			slurmFile.println();

			lsfFile.println("#!/usr/bin/tcsh");
			lsfFile.println("#BSUB -P "+acc.trim());
			lsfFile.println("#BSUB -J "+jobName.trim());
			lsfFile.println("#BSUB -W "+limit.trim());
			lsfFile.println("#BSUB -q batch ");
			lsfFile.println("#BSUB -o job.o%J");
			lsfFile.println("#BSUB -e job.e%J");
			lsfFile.println("#BSUB -alloc_flags \"gpumps\"");
			lsfFile.println("#BSUB -nnodes "+ranks);
			lsfFile.println();
			lsfFile.println("#module purge");
			lsfFile.println("#module load intel openmpi");
			lsfFile.println();
			lsfFile.println("#source required files");
			lsfFile.println();
			lsfFile.println("#cd to job directory");
			lsfFile.println();

			runFile.println("# OpenMolcas runs for the spin states are run in the same workspace "+fp.trim()+fs.trim());
			slurmFile.println("# OpenMolcas runs for the spin states are run in the same workspace "+fp.trim()+fs.trim());
			lsfFile.println("# OpenMolcas runs for the spin states are run in the same workspace "+fp.trim()+fs.trim());

			runFile.println("# Evaluate the integrals");
			slurmFile.println("# Evaluate the integrals");
			lsfFile.println("# Evaluate the integrals");
			fullName = fp.trim()+fs.trim()+"_INT";
			runFile.println("cp "+fullName.trim()+".input "+fp.trim()+fs.trim()+".input; "+"pymolcas "+fp.trim()+fs.trim()+".input > "+fullName.trim()+".output");
			slurmFile.println("cp "+fullName.trim()+".input "+fp.trim()+fs.trim()+".input; "+"pymolcas "+fp.trim()+fs.trim()+".input > "+fullName.trim()+".output");
			lsfFile.println("cp "+fullName.trim()+".input "+fp.trim()+fs.trim()+".input; "+"pymolcas "+fp.trim()+fs.trim()+".input > "+fullName.trim()+".output");

			runFile.println("# Remove the generic input file to avoid clutter");
			slurmFile.println("# Remove the generic input file to avoid clutter");
			lsfFile.println("# Remove the generic input file to avoid clutter");
			runFile.println("rm "+fp.trim()+fs.trim()+".input");
			slurmFile.println("rm "+fp.trim()+fs.trim()+".input");
			lsfFile.println("rm "+fp.trim()+fs.trim()+".input");
			
			runFile.println("# OpenMolcas restricts rhe runfilename to 8 characters requiring copying it to the generic RUNFILE");
			slurmFile.println("# OpenMolcas restricts rhe runfilename to 8 characters requiring copying it to the generic RUNFILE");
			lsfFile.println("# OpenMolcas restricts rhe runfilename to 8 characters requiring copying it to the generic RUNFILE");
			runFile.println("cp "+fp.trim()+fs.trim()+".runfil RUNFILE");
			slurmFile.println("cp "+fp.trim()+fs.trim()+".runfil RUNFILE");
			lsfFile.println("cp "+fp.trim()+fs.trim()+".runfil RUNFILE");
			runFile.println("# Rotation of all spin states");
			slurmFile.println("# Rotation of all spin states");
			lsfFile.println("# Rotation of all spin states");
			runFile.println("grotate < "+fp.trim()+name1+"_rotate.input");
//			for(int i=0; i<numStates; i++) {
//				Integer index = ndxStateList[stateset][i];
//				runFile.println("cp "+fp.trim()+name2+"_"+stateNames[index]+".INPORB INPORB");
//				runFile.println("rotharm < "+fp.trim()+name1+"_rotate.input");
//				runFile.println("grotate < "+fp.trim()+name1+"_rotate.input");
//				runFile.println("mv ROTORB "+fp.trim()+name1+"_"+stateNames[index]+".ROTORB");
//				runFile.println("mv transrot.xyz "+fp.trim()+name1+".xyz");
//				runFile.println("cp "+fp.trim()+name2+"_"+stateNames[index]+".det "+fp.trim()+name1+"_"+stateNames[index]+".det");
//				slurmFile.println("cp "+fp.trim()+name2+"_"+stateNames[index]+".INPORB INPORB");
//				slurmFile.println("rotharm < "+fp.trim()+name1+"_rotate.input");
//				slurmFile.println("grotate < "+fp.trim()+name1+"_rotate.input");
//				slurmFile.println("mv ROTORB "+fp.trim()+name1+"_"+stateNames[index]+".ROTORB");
//				slurmFile.println("mv transrot.xyz "+fp.trim()+name1+".xyz");
//				slurmFile.println("cp "+fp.trim()+name2+"_"+stateNames[index]+".det "+fp.trim()+name1+"_"+stateNames[index]+".det");
//				lsfFile.println("cp "+fp.trim()+name2+"_"+stateNames[index]+".INPORB INPORB");
//				lsfFile.println("rotharm < "+fp.trim()+name1+"_rotate.input");
//				lsfFile.println("grotate < "+fp.trim()+name1+"_rotate.input");
//				lsfFile.println("mv ROTORB "+fp.trim()+name1+"_"+stateNames[index]+".ROTORB");
//				lsfFile.println("mv transrot.xyz "+fp.trim()+name1+".xyz");
//				lsfFile.println("cp "+fp.trim()+name2+"_"+stateNames[index]+".det "+fp.trim()+name1+"_"+stateNames[index]+".det");
//			}
			runFile.println("# Copy the determinant file of reference fragment as it remains unchanged upon rotation");
			slurmFile.println("# Copy the determinant file of reference fragment as it remains unchanged upon rotation");
			lsfFile.println("# Copy the determinant file of reference fragment as it remains unchanged upon rotation");
			for(int i=0; i<numStates; i++) {
				Integer index = ndxStateList[stateset][i];
//				runFile.println("cp "+fp.trim()+name2+"_"+stateNames[index]+".INPORB INPORB");
//				runFile.println("rotharm < "+fp.trim()+name1+"_rotate.input");
//				runFile.println("grotate < "+fp.trim()+name1+"_rotate.input");
//				runFile.println("mv ROTORB "+fp.trim()+name1+"_"+stateNames[index]+".ROTORB");
//				runFile.println("mv transrot.xyz "+fp.trim()+name1+".xyz");
				runFile.println("cp "+fp.trim()+name2+"_"+stateNames[index]+".det "+fp.trim()+name1+"_"+stateNames[index]+".det");
//				slurmFile.println("cp "+fp.trim()+name2+"_"+stateNames[index]+".INPORB INPORB");
//				slurmFile.println("rotharm < "+fp.trim()+name1+"_rotate.input");
//				slurmFile.println("grotate < "+fp.trim()+name1+"_rotate.input");
//				slurmFile.println("mv ROTORB "+fp.trim()+name1+"_"+stateNames[index]+".ROTORB");
//				slurmFile.println("mv transrot.xyz "+fp.trim()+name1+".xyz");
				slurmFile.println("cp "+fp.trim()+name2+"_"+stateNames[index]+".det "+fp.trim()+name1+"_"+stateNames[index]+".det");
//				lsfFile.println("cp "+fp.trim()+name2+"_"+stateNames[index]+".INPORB INPORB");
//				lsfFile.println("rotharm < "+fp.trim()+name1+"_rotate.input");
//				lsfFile.println("grotate < "+fp.trim()+name1+"_rotate.input");
//				lsfFile.println("mv ROTORB "+fp.trim()+name1+"_"+stateNames[index]+".ROTORB");
//				lsfFile.println("mv transrot.xyz "+fp.trim()+name1+".xyz");
				lsfFile.println("cp "+fp.trim()+name2+"_"+stateNames[index]+".det "+fp.trim()+name1+"_"+stateNames[index]+".det");
			}
		runFile.close();
		slurmFile.close();
		lsfFile.close();
		} catch(IOException ei) {
		}
	}
	
	public void create_Basis(String basis, String contract) {
		for(int i=0; i<elements.length; i++) {
			if(elements[i].equals("H")) {
				if(basis.equals("ano-s") && contract.equals("s")) basisSet[i]="h.ano-s...2s1p.";
				if(basis.equals("ano-s") && contract.equals("m")) basisSet[i]="h.ano-s...3s2p.";
				if(basis.equals("ano-s") && contract.equals("l")) basisSet[i]="h.ano-s...4s3p.";
				if(basis.equals("ano-l") && contract.equals("s")) basisSet[i]="h.ano-l...3s2p.";
				if(basis.equals("ano-l") && contract.equals("m")) basisSet[i]="h.ano-l...3s2p1d.";
				if(basis.equals("ano-l") && contract.equals("l")) basisSet[i]="h.ano-l...6s4p3d.";
			}
			if(elements[i].equals("C")) {
				if(basis.equals("ano-s") && contract.equals("s")) basisSet[i]="c.ano-s...3s2p1d.";
				if(basis.equals("ano-s") && contract.equals("m")) basisSet[i]="c.ano-s...4s3p2d.";
				if(basis.equals("ano-s") && contract.equals("l")) basisSet[i]="c.ano-s...7s6p3d.";
				if(basis.equals("ano-l") && contract.equals("s")) basisSet[i]="c.ano-l...4s3p2d.";
				if(basis.equals("ano-l") && contract.equals("m")) basisSet[i]="c.ano-l...5s4p3d2f.";
				if(basis.equals("ano-l") && contract.equals("l")) basisSet[i]="c.ano-l...7s6p4d3f.";
			}
			if(elements[i].equals("N")) {
				if(basis.equals("ano-s") && contract.equals("s")) basisSet[i]="n.ano-s...3s2p1d.";
				if(basis.equals("ano-s") && contract.equals("m")) basisSet[i]="n.ano-s...4s3p2d.";
				if(basis.equals("ano-s") && contract.equals("l")) basisSet[i]="n.ano-s...7s6p3d.";
				if(basis.equals("ano-l") && contract.equals("s")) basisSet[i]="n.ano-l...4s3p2d.";
				if(basis.equals("ano-l") && contract.equals("m")) basisSet[i]="n.ano-l...5s4p3d2f.";
				if(basis.equals("ano-l") && contract.equals("l")) basisSet[i]="n.ano-l...7s6p4d3f.";
			}
			if(elements[i].equals("O")) {
				if(basis.equals("ano-s") && contract.equals("s")) basisSet[i]="o.ano-s...3s2p1d.";
				if(basis.equals("ano-s") && contract.equals("m")) basisSet[i]="o.ano-s...4s3p2d.";
				if(basis.equals("ano-s") && contract.equals("l")) basisSet[i]="o.ano-s...7s6p3d.";
				if(basis.equals("ano-l") && contract.equals("s")) basisSet[i]="o.ano-l...4s3p2d.";
				if(basis.equals("ano-l") && contract.equals("m")) basisSet[i]="o.ano-l...5s4p3d2f.";
				if(basis.equals("ano-l") && contract.equals("l")) basisSet[i]="o.ano-l...7s6p4d3f.";
			}
			if(elements[i].equals("F")) {
				if(basis.equals("ano-s") && contract.equals("s")) basisSet[i]="f.ano-s...3s2p1d.";
				if(basis.equals("ano-s") && contract.equals("m")) basisSet[i]="f.ano-s...4s3p2d.";
				if(basis.equals("ano-s") && contract.equals("l")) basisSet[i]="f.ano-s...7s6p3d.";
				if(basis.equals("ano-l") && contract.equals("s")) basisSet[i]="f.ano-l...4s3p2d.";
				if(basis.equals("ano-l") && contract.equals("m")) basisSet[i]="f.ano-l...5s4p3d2f.";
				if(basis.equals("ano-l") && contract.equals("l")) basisSet[i]="f.ano-l...7s6p4d3f.";
			}
			if(elements[i].equals("P")) {
				if(basis.equals("ano-s") && contract.equals("s")) basisSet[i]="p.ano-s...4s3p2d.";
				if(basis.equals("ano-s") && contract.equals("m")) basisSet[i]="p.ano-s...5s4p3d.";
				if(basis.equals("ano-s") && contract.equals("l")) basisSet[i]="p.ano-s...7s7p4d.";
				if(basis.equals("ano-l") && contract.equals("s")) basisSet[i]="p.ano-l...5s4p2d.";
				if(basis.equals("ano-l") && contract.equals("m")) basisSet[i]="p.ano-l...6s5p4d3f.";
				if(basis.equals("ano-l") && contract.equals("l")) basisSet[i]="p.ano-l...7s7p5d4f.";
			}
			if(elements[i].equals("S")) {
				if(basis.equals("ano-s") && contract.equals("s")) basisSet[i]="s.ano-s...4s3p2d.";
				if(basis.equals("ano-s") && contract.equals("m")) basisSet[i]="s.ano-s...5s4p3d.";
				if(basis.equals("ano-s") && contract.equals("l")) basisSet[i]="s.ano-s...7s7p4d.";
				if(basis.equals("ano-l") && contract.equals("s")) basisSet[i]="s.ano-l...5s4p2d.";
				if(basis.equals("ano-l") && contract.equals("m")) basisSet[i]="s.ano-l...6s5p4d3f.";
				if(basis.equals("ano-l") && contract.equals("l")) basisSet[i]="s.ano-l...7s7p5d4f.";
			}
			if(elements[i].equals("Cl")) {
				if(basis.equals("ano-s") && contract.equals("s")) basisSet[i]="cl.ano-s...4s3p2d.";
				if(basis.equals("ano-s") && contract.equals("m")) basisSet[i]="cl.ano-s...5s4p3d.";
				if(basis.equals("ano-s") && contract.equals("l")) basisSet[i]="cl.ano-s...7s7p4d.";
				if(basis.equals("ano-l") && contract.equals("s")) basisSet[i]="cl.ano-l...5s4p2d.";
				if(basis.equals("ano-l") && contract.equals("m")) basisSet[i]="cl.ano-l...6s5p4d3f.";
				if(basis.equals("ano-l") && contract.equals("l")) basisSet[i]="cl.ano-l...7s7p5d4f.";
			}
			if(elements[i].equals("Se")) {
				if(basis.equals("ano-s") && contract.equals("s")) basisSet[i]="se.ano-s...5s4p3d.";
				if(basis.equals("ano-s") && contract.equals("m")) basisSet[i]="se.ano-s...6s5p4d.";
				if(basis.equals("ano-s") && contract.equals("l")) basisSet[i]="se.ano-s...9s9p5d.";
			}
			if(elements[i].equals("Br")) {
				if(basis.equals("ano-s") && contract.equals("s")) basisSet[i]="br.ano-s...5s4p3d.";
				if(basis.equals("ano-s") && contract.equals("m")) basisSet[i]="br.ano-s...6s5p4d.";
				if(basis.equals("ano-s") && contract.equals("l")) basisSet[i]="br.ano-s...9s9p5d.";
			}
		}
	}
	
	public Integer basis_Index(String atom) {
		Integer ndx=0;
		for(int i=0; i<elements.length; i++) {
			if(atom.trim().equals(elements[i].trim())) ndx=i;
		}
		return ndx;
	}
	
	public void write_Run_Script_MEBFs(String p, String pn, Integer nfrags, String[] fname, String[] frags, String[] source, Integer[] fstat, Integer[] lenStateList, Integer[][] ndxStateList, Integer ranks, Integer memory, Object[][] fragDef, String acc, String jobName, String limit) {
		String fileName = p+"_Molcas.run";
		String slurmName = p+"_Molcas.slurm";
		String lsfName = p+"_Molcas.lsf";
		String fullName;
		String fullName2;
		try {
			PrintfWriter runFile = new PrintfWriter(new FileWriter(fileName));
			PrintfWriter slurmFile = new PrintfWriter(new FileWriter(slurmName));
			PrintfWriter lsfFile = new PrintfWriter(new FileWriter(lsfName));
			runFile.println("#!/usr/bin/tcsh");
			slurmFile.println("#!/usr/bin/tcsh");
			lsfFile.println("#!/usr/bin/tcsh");
			
			runFile.println("setenv MOLCAS_NPROCS "+ranks);
			runFile.println("setenv MOLCAS_MEM "+memory);
			slurmFile.println("setenv MOLCAS_NPROCS "+ranks);
			slurmFile.println("setenv MOLCAS_MEM "+memory);
			lsfFile.println("setenv MOLCAS_NPROCS "+ranks);
			lsfFile.println("setenv MOLCAS_MEM "+memory);

			slurmFile.println("#SBATCH --account="+acc.trim());
			slurmFile.println("#SBATCH --nodes=1");
			slurmFile.println("#SBATCH --tasks-per-node="+ranks);
			slurmFile.println("#SBATCH --job-name="+jobName.trim());
			slurmFile.println("#SBATCH --time="+limit.trim());
			slurmFile.println();
			slurmFile.println("module purge");
			slurmFile.println("module load intel openmpi");
			slurmFile.println();
			slurmFile.println("#source required files");
			slurmFile.println();
			slurmFile.println("#cd to job directory");
			slurmFile.println();
			slurmFile.println("setenv MOLCAS_NPROCS "+ranks);
			slurmFile.println("setenv MOLCAS_MEM "+memory);
			
			lsfFile.println("#BSUB -P "+acc.trim());
			lsfFile.println("#BSUB -J "+jobName.trim());
			lsfFile.println("#BSUB -W "+limit.trim());
			lsfFile.println("#BSUB -q batch ");
			lsfFile.println("#BSUB -o job.o%J");
			lsfFile.println("#BSUB -e job.e%J");
			lsfFile.println("#BSUB -alloc_flags \"gpumps\"");
			lsfFile.println("#BSUB -nnodes "+ranks);
			lsfFile.println();
			lsfFile.println("#module purge");
			lsfFile.println("#module load intel openmpi");
			lsfFile.println();
			lsfFile.println("#source required files");
			lsfFile.println();
			lsfFile.println("#cd to job directory");
			lsfFile.println();
			lsfFile.println("setenv MOLCAS_NPROCS "+ranks);
			lsfFile.println("setenv MOLCAS_MEM "+memory);
			
			runFile.println("# Merge the coordinates of the fragments into a single XYZ file "+p.trim()+".xyz");
			slurmFile.println("# Merge the coordinates of the fragments into a single XYZ file "+p.trim()+".xyz");
			lsfFile.println("# Merge the coordinates of the fragments into a single XYZ file "+p.trim()+".xyz");
			
			runFile.print("merge_xyz "+p.trim());
			slurmFile.print("merge_xyz "+p.trim());
			lsfFile.print("merge_xyz "+p.trim());
			for(int i=0; i<nfrags; i++) {
				runFile.print(" "+pn.trim()+frags[i].trim());
				slurmFile.print(" "+pn.trim()+frags[i].trim());
				lsfFile.print(" "+pn.trim()+frags[i].trim());
			}
			runFile.println();
			slurmFile.println();
			lsfFile.println();
			runFile.println("# One electron integrals for this MEBF in the workspace "+p.trim());
			slurmFile.println("# One electron integrals for this MEBF in the workspace "+p.trim());
			lsfFile.println("# One electron integrals for this MEBF in the workspace "+p.trim());
			fullName = p.trim()+"_ONE";
			runFile.println("cp "+fullName.trim()+".input "+p.trim()+".input; "+"pymolcas "+p.trim()+".input > "+fullName.trim()+".output");
			slurmFile.println("cp "+fullName.trim()+".input "+p.trim()+".input; "+"pymolcas "+p.trim()+".input > "+fullName.trim()+".output");
			lsfFile.println("cp "+fullName.trim()+".input "+p.trim()+".input; "+"pymolcas "+p.trim()+".input > "+fullName.trim()+".output");

			runFile.println("# Remove the generic input file to avoid clutter");
			slurmFile.println("# Remove the generic input file to avoid clutter");
			lsfFile.println("# Remove the generic input file to avoid clutter");
			runFile.println("rm "+p.trim()+".input");
			slurmFile.println("rm "+p.trim()+".input");
			lsfFile.println("rm "+p.trim()+".input");
						
			Integer index=0;
			Integer stateIndex = 0;
			for(int i=0; i<nfrags; i++) {
				Boolean rot=!frags[i].equals(source[i]);
				runFile.println("cp "+pn.trim()+frags[i].trim()+".runfil RUNFIL"+frags[i]);
				runFile.println("cp "+pn.trim()+frags[i].trim()+".oneint ONEINT"+frags[i]);
				slurmFile.println("cp "+pn.trim()+frags[i].trim()+".runfil RUNFIL"+frags[i]);
				slurmFile.println("cp "+pn.trim()+frags[i].trim()+".oneint ONEINT"+frags[i]);
				lsfFile.println("cp "+pn.trim()+frags[i].trim()+".runfil RUNFIL"+frags[i]);
				lsfFile.println("cp "+pn.trim()+frags[i].trim()+".oneint ONEINT"+frags[i]);
//		for(int j=0; j<lenStateList[fstat[i]]; j++) {
//			index++;
//			stateIndex=ndxStateList[fstat[i]][j];
//			if(rot) {
//				runFile.println("cp "+pn.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".ROTORB INPORB."+(i+1)+"_"+(j+1));
//				slurmFile.println("cp "+pn.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".ROTORB INPORB."+(i+1)+"_"+(j+1));
//				lsfFile.println("cp "+pn.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".ROTORB INPORB."+(i+1)+"_"+(j+1));
//			} else {
//				runFile.println("cp "+pn.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".INPORB INPORB."+(i+1)+"_"+(j+1));
//				slurmFile.println("cp "+pn.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".INPORB INPORB."+(i+1)+"_"+(j+1));
//				lsfFile.println("cp "+pn.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".INPORB INPORB."+(i+1)+"_"+(j+1));
//			};
//			String sIndex="00"+index;
//			if(index>9) sIndex="0"+index;
//			runFile.println("cp "+pn.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".det "+p.trim()+"_"+sIndex.trim()+".det");
//			slurmFile.println("cp "+pn.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".det "+p.trim()+"_"+sIndex.trim()+".det");
//			lsfFile.println("cp "+pn.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".det "+p.trim()+"_"+sIndex.trim()+".det");
//					if(index<10) {
//						runFile.println("cp "+pn.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".det "+p.trim()+"_00"+index+".det");
//					} else {
//						runFile.println("cp "+pn.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".det "+p.trim()+"_0"+index+".det");
//					}
//		}
			}
			runFile.println("cp "+p.trim()+".runfil RUNFIL");
			slurmFile.println("cp "+p.trim()+".runfil RUNFIL");
			lsfFile.println("cp "+p.trim()+".runfil RUNFIL");
			fullName = p.trim()+"_CB";
//			runFile.println("cp "+fullName.trim()+".input "+p.trim()+".input; "+"common_basis < "+p.trim()+".input > "+fullName.trim()+".output");
//			slurmFile.println("cp "+fullName.trim()+".input "+p.trim()+".input; "+"common_basis < "+p.trim()+".input > "+fullName.trim()+".output");
//			lsfFile.println("cp "+fullName.trim()+".input "+p.trim()+".input; "+"common_basis < "+p.trim()+".input > "+fullName.trim()+".output");
			runFile.println("gcommon < "+fullName.trim()+".input > "+fullName.trim()+".output");
			slurmFile.println("gcommon < "+fullName.trim()+".input > "+fullName.trim()+".output");
			lsfFile.println("gcommon < "+fullName.trim()+".input > "+fullName.trim()+".output");
			index=0;
			for(int i=0; i<nfrags; i++) {
				runFile.println("rm RUNFIL"+frags[i]);
				slurmFile.println("rm RUNFIL"+frags[i]);
				lsfFile.println("rm RUNFIL"+frags[i]);
			}
			
//				runFile.println("rm ONEINT"+(i+1));
//				for(int j=0; j<lenStateList[fstat[i]]; j++) {
//					index++;
//					stateIndex=ndxStateList[fstat[i]][j];
//					if(index<10) {
//						runFile.println("mv "+p.trim()+"_00"+index+".vec "+p.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".vec");
//						runFile.println("mv "+p.trim()+"_00"+index+".det "+p.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".det");
//						slurmFile.println("mv "+p.trim()+"_00"+index+".vec "+p.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".vec");
//						slurmFile.println("mv "+p.trim()+"_00"+index+".det "+p.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".det");
//						lsfFile.println("mv "+p.trim()+"_00"+index+".vec "+p.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".vec");
//						lsfFile.println("mv "+p.trim()+"_00"+index+".det "+p.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".det");
//					} else {
//						runFile.println("mv "+p.trim()+"_0"+index+".vec "+p.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".vec");
//						runFile.println("mv "+p.trim()+"_0"+index+".det "+p.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".det");
//						slurmFile.println("mv "+p.trim()+"_0"+index+".vec "+p.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".vec");
//						slurmFile.println("mv "+p.trim()+"_0"+index+".det "+p.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".det");
//						lsfFile.println("mv "+p.trim()+"_0"+index+".vec "+p.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".vec");
//						lsfFile.println("mv "+p.trim()+"_0"+index+".det "+p.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".det");
//					}
//					runFile.println("det_header "+pn.trim()+source[i].trim()+"_"+stateNames[stateIndex].trim()+".output "+p.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".det");
//					runFile.println("rm INPORB."+(i+1)+"_"+(j+1));
//					slurmFile.println("det_header "+pn.trim()+source[i].trim()+"_"+stateNames[stateIndex].trim()+".output "+p.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".det");
//					slurmFile.println("rm INPORB."+(i+1)+"_"+(j+1));
//					lsfFile.println("det_header "+pn.trim()+source[i].trim()+"_"+stateNames[stateIndex].trim()+".output "+p.trim()+frags[i].trim()+"_"+stateNames[stateIndex].trim()+".det");
//					lsfFile.println("rm INPORB."+(i+1)+"_"+(j+1));
//				}	
//			}
			
			fullName = p.trim()+"_CB";
			runFile.println("setenv DELETED ` grep \"Deleted orbitals in MOTRA\" "+fullName.trim()+".output | cut -b 30-34 `");
			runFile.println("touch TRAINT");
			slurmFile.println("setenv DELETED ` grep \"Deleted orbitals in MOTRA\" "+fullName.trim()+".output | cut -b 30-34 `");
			slurmFile.println("touch TRAINT");
			lsfFile.println("setenv DELETED ` grep \"Deleted orbitals in MOTRA\" "+fullName.trim()+".output | cut -b 30-34 `");
			lsfFile.println("touch TRAINT");
			fullName = p.trim()+"_TWO";
			runFile.println("cp "+fullName.trim()+".input "+p.trim()+".input; "+"pymolcas "+p.trim()+".input > "+fullName.trim()+".output");
			runFile.println("setenv OMP_NUM_THREADS 12");
			runFile.println("cp "+p.trim()+".RUNFILE RUNFILE");
			runFile.println("cp "+p.trim()+"_CHMOT1 _CHMOT1");
			runFile.println("cp "+p.trim()+".CHORST CHORST");
			runFile.println("cp "+p.trim()+".ONEINT ONEINT");
			runFile.println("cp "+p.trim()+".TRAONE TRAONE");
			runFile.println("cp "+p.trim()+".CHOMAP CHOMAP");
			runFile.println("cp "+p.trim()+".CHRED  CHRED");
			runFile.println("cp "+p.trim()+".CHVEC1 CHVEC1");
			runFile.println("cp "+p.trim()+".comorb COMMONORB");
			slurmFile.println("cp "+fullName.trim()+".input "+p.trim()+".input; "+"pymolcas "+p.trim()+".input > "+fullName.trim()+".output");
			slurmFile.println("setenv OMP_NUM_THREADS "+ranks);
			slurmFile.println("cp "+p.trim()+".RUNFILE RUNFILE");
			slurmFile.println("cp "+p.trim()+"_CHMOT1 <");
			slurmFile.println("cp "+p.trim()+".CHORST CHORST");
			slurmFile.println("cp "+p.trim()+".ONEINT ONEINT");
			slurmFile.println("cp "+p.trim()+".TRAONE TRAONE");
			slurmFile.println("cp "+p.trim()+".CHOMAP CHOMAP");
			slurmFile.println("cp "+p.trim()+".CHRED  CHRED");
			slurmFile.println("cp "+p.trim()+".CHVEC1 CHVEC1");
			slurmFile.println("cp "+p.trim()+".comorb COMMONORB");
			lsfFile.println("cp "+fullName.trim()+".input "+p.trim()+".input; "+"pymolcas "+p.trim()+".input > "+fullName.trim()+".output");
			lsfFile.println("setenv OMP_NUM_THREADS "+ranks);
//			lsfFile.println("cp "+p.trim()+".runfil RUNFILE");
			lsfFile.println("cp "+p.trim()+"_CHMOT1 _CHMOT1");
			lsfFile.println("cp "+p.trim()+".CHORST CHORST");
			lsfFile.println("cp "+p.trim()+".ONEINT ONEINT");
			lsfFile.println("cp "+p.trim()+".TRAONE TRAONE");
			lsfFile.println("cp "+p.trim()+".CHOMAP CHOMAP");
			lsfFile.println("cp "+p.trim()+".CHRED  CHRED");
			lsfFile.println("cp "+p.trim()+".CHVEC1 CHVEC1");
			lsfFile.println("cp "+p.trim()+".comorb COMMONORB");
			fullName = p.trim()+"_RDCH";
			runFile.println("rdcho $MOLCAS_NPROCS > "+fullName.trim()+".output");
			runFile.println("rm _CHMOT1");
//			runFile.println("cp "+p.trim()+".RUNFILE RUNFILE");
//			runFile.println("cp "+p.trim()+".CHORST CHORST");
//			runFile.println("cp "+p.trim()+".ONEINT ONEINT");
//			runFile.println("cp "+p.trim()+".TRAONE TRAONE");
//			runFile.println("cp "+p.trim()+".CHOMAP CHOMAP");
//			runFile.println("cp "+p.trim()+".CHRED  CHRED");
//			runFile.println("cp "+p.trim()+".CHVEC1 CHVEC1");
//			runFile.println("cp "+p.trim()+".comorb COMMONORB");
			slurmFile.println("rdcho $MOLCAS_NPROCS > "+fullName.trim()+".output");
			slurmFile.println("rm _CHMOT1");
//			slurmFile.println("cp "+p.trim()+".RUNFILE RUNFILE");
//			slurmFile.println("cp "+p.trim()+".CHORST CHORST");
//			slurmFile.println("cp "+p.trim()+".ONEINT ONEINT");
//			slurmFile.println("cp "+p.trim()+".TRAONE TRAONE");
//			slurmFile.println("cp "+p.trim()+".CHOMAP CHOMAP");
//			slurmFile.println("cp "+p.trim()+".CHRED  CHRED");
//			slurmFile.println("cp "+p.trim()+".CHVEC1 CHVEC1");
//			slurmFile.println("cp "+p.trim()+".comorb COMMONORB");
			lsfFile.println("rdcho $MOLCAS_NPROCS > "+fullName.trim()+".output");
			lsfFile.println("rm _CHMOT1");
//			lsfFile.println("cp "+p.trim()+".RUNFILE RUNFILE");
//			lsfFile.println("cp "+p.trim()+".CHORST CHORST");
//			lsfFile.println("cp "+p.trim()+".ONEINT ONEINT");
//			lsfFile.println("cp "+p.trim()+".TRAONE TRAONE");
//			lsfFile.println("cp "+p.trim()+".CHOMAP CHOMAP");
//			lsfFile.println("cp "+p.trim()+".CHRED  CHRED");
//			lsfFile.println("cp "+p.trim()+".CHVEC1 CHVEC1");
//			lsfFile.println("cp "+p.trim()+".comorb COMMONORB");
			fullName = p.trim()+"_CB";
			fullName2 = p.trim()+"_MEBFRDTR";
			runFile.println("rdtraint < "+fullName.trim()+".input > "+fullName2+".output");
			runFile.println("mv COMMONORB "+p.trim()+".comorb");
			runFile.println("rm RUNFILE");
			runFile.println("rm CHORST");
			runFile.println("rm ONEINT*");
			runFile.println("rm TRAONE");
			runFile.println("rm TRAINT");
			runFile.println("rm CHOMAP");
			runFile.println("rm CHRED");
			runFile.println("rm CHVEC1");
			runFile.println("unsetenv DELETED");
			slurmFile.println("rdtraint < "+fullName.trim()+".input > "+fullName2+".output");
			slurmFile.println("mv COMMONORB "+p.trim()+".comorb");
			slurmFile.println("rm RUNFILE");
			slurmFile.println("rm CHORST");
			slurmFile.println("rm ONEINT*");
			slurmFile.println("rm TRAONE");
			slurmFile.println("rm TRAINT");
			slurmFile.println("rm CHOMAP");
			slurmFile.println("rm CHRED");
			slurmFile.println("rm CHVEC1");
			slurmFile.println("unsetenv DELETED");
			lsfFile.println("rdtraint < "+fullName.trim()+".input > "+fullName2+".output");
			lsfFile.println("mv COMMONORB "+p.trim()+".comorb");
			lsfFile.println("rm RUNFILE");
			lsfFile.println("rm CHORST");
			lsfFile.println("rm ONEINT*");
			lsfFile.println("rm TRAONE");
			lsfFile.println("rm TRAINT");
			lsfFile.println("rm CHOMAP");
			lsfFile.println("rm CHRED");
			lsfFile.println("rm CHVEC1");
			lsfFile.println("unsetenv DELETED");
			runFile.close();
			slurmFile.close();
			lsfFile.close();
		} catch(IOException ei) {
		}
		fileName = p+"_GronOR.run";
		slurmName = p+"_GronOR.slurm";
		lsfName = p+"_GronOR.lsf";
		try {
			PrintfWriter runFile = new PrintfWriter(new FileWriter(fileName));
			PrintfWriter slurmFile = new PrintfWriter(new FileWriter(slurmName));
			PrintfWriter lsfFile = new PrintfWriter(new FileWriter(lsfName));
			runFile.println("#!/usr/bin/tcsh");
			slurmFile.println("#!/usr/bin/tcsh");
			lsfFile.println("#!/usr/bin/tcsh");

			slurmFile.println("#SBATCH --account="+acc.trim());
			slurmFile.println("#SBATCH --nodes=1");
			slurmFile.println("#SBATCH --tasks-per-node="+ranks);
			slurmFile.println("#SBATCH --job-name="+jobName.trim());
			slurmFile.println("#SBATCH --time="+limit.trim());
			slurmFile.println();
			slurmFile.println("module purge");
			slurmFile.println("module load intel openmpi");
			slurmFile.println();
			slurmFile.println("#source required files");
			slurmFile.println();
			slurmFile.println("#cd to job directory");
			slurmFile.println();
			
			lsfFile.println("#BSUB -P "+acc.trim());
			lsfFile.println("#BSUB -J "+jobName.trim());
			lsfFile.println("#BSUB -W "+limit.trim());
			lsfFile.println("#BSUB -q batch ");
			lsfFile.println("#BSUB -o job.o%J");
			lsfFile.println("#BSUB -e job.e%J");
			lsfFile.println("#BSUB -alloc_flags \"gpumps\"");
			lsfFile.println("#BSUB -nnodes "+ranks);
			lsfFile.println();
			lsfFile.println("#module purge");
			lsfFile.println("#module load ");
			lsfFile.println();
			lsfFile.println("#source required files");
			lsfFile.println();
			lsfFile.println("#cd to job directory");
			lsfFile.println();
			fullName= p+"_GronOR";
			runFile.println("mpirun -n "+ranks+" gronor "+fullName);
			slurmFile.println("mpirun -n "+ranks+" gronor "+fullName);
			lsfFile.println("cp /gpfs/alpine/chm154/proj-shared/gronor .");
			lsfFile.println();
			lsfFile.println("jsrun -a1 -g1 -r6 -c7 -d cyclic ./gronor "+fullName);
			runFile.close();
			slurmFile.close();
			lsfFile.close();
		} catch(IOException ei) {
		}
	}
	
	public Boolean read_XYZ() {
		String fileName = fragmentName+".xyz";
		if(fragmentName.trim().equals("")) return false;
		String card;
		StringTokenizer st;
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			card=br.readLine();
			st = new StringTokenizer(card," ");
			numAtoms=Integer.valueOf(st.nextToken());
//			System.out.println("READING "+fileName+" : "+numAtoms);
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
//			System.out.println("Read:    "+fileName);
			return true;
		} catch(IOException ef) {
			return false;
		}
	}

	public Integer read_Alt(String nameP, String nameA) {
		String fileName = nameP.trim()+nameA.trim()+".alter";
		String card;
		numAlt=0;
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			card=br.readLine();
			numAlt=Integer.valueOf(card.substring(0,6).trim());
			
			for(int i=0; i<numAlt; i++) {
				card=br.readLine();
				alter[i][0]=Integer.valueOf(card.substring(0,6).trim());
				alter[i][1]=Integer.valueOf(card.substring(6,12).trim());
			}
			br.close();
			return numAlt;
		} catch(IOException ef) {
			return numAlt;
		}
	}
	
	public Boolean write_XYZ() {
		String fileName = fragmentName+".xyz";
//		System.out.println("WRITING "+fileName+" : "+numAtoms);
		File f = new File(fileName);
		Boolean skip=true;
		for(int i=0; i<6; i++) if(RandT[i]!=0.0) skip=false;
		if(!f.exists()) skip=false;
		if(!skip) {
			try {
				PrintfWriter xyzFile = new PrintfWriter(new FileWriter(fileName));
				xyzFile.printf("%6d",numAtoms);
				xyzFile.println();
				xyzFile.println("Coordinates in Angstrom");
				for(int i=0; i<numAtoms; i++) {
					xyzFile.print(atomLabel[i]); 
					xyzFile.printf("%16.8f",coordinates[i][0]);
					xyzFile.printf("%16.8f",coordinates[i][1]);
					xyzFile.printf("%16.8f",coordinates[i][2]);
					xyzFile.println();
				}
				xyzFile.close();
//				System.out.println("Written: "+fileName);
//				return true;
			} catch(IOException ei) {
				return false;
			}
		}

		if(convergedName.trim().length()>1) {
			try {
				PrintfWriter xyzFile = new PrintfWriter(new FileWriter(convergedName));
				convergedName=" ";
				xyzFile.printf("%6d",numAtoms);
				xyzFile.println();
				xyzFile.println("DFT converged coordinates in Angstrom");
				for(int i=0; i<numAtoms; i++) {
					xyzFile.print(atomLabel[i]); 
					xyzFile.printf("%16.8f",coordinates[i][0]);
					xyzFile.printf("%16.8f",coordinates[i][1]);
					xyzFile.printf("%16.8f",coordinates[i][2]);
					xyzFile.println();
				}
				xyzFile.close();
//				return true;
			} catch(IOException ei) {
				return false;
			}
		}
		
		return true;
	}
	
	public Boolean write_rotharm_input(String fp, String na, String nb, Integer stateset, Integer[] lenStateList, Integer[][] ndxStateList) {
		String fileName = fp.trim()+na.trim()+"_rotate.input";
		Integer numStates = lenStateList[stateset];
//		System.out.println(fileName);
		/*
		double dmax;
		dmax=0.0;
		int imax=0, jmax=0, kmax=0;
		for(int i=0; i<numAtoms-2; i++) {
			if(!atomLabel[i].trim().equals("H")) {
				for(int j=i+1; j<numAtoms-1; j++) {
					if(!atomLabel[j].trim().equals("H")) {
						double dij=0.0;
						for(int m=0; m<3; m++) dij=dij+(coordinates[i][m]-coordinates[j][m])*(coordinates[i][m]-coordinates[j][m]);
						for(int k=j+1; k<numAtoms; k++) {
							if(!atomLabel[k].trim().equals("H")) {
								double dik=0.0;
								double djk=0.0;
								for(int m=0; m<3; m++) dik=dik+(coordinates[i][m]-coordinates[k][m])*(coordinates[i][m]-coordinates[k][m]);
								for(int m=0; m<3; m++) djk=djk+(coordinates[j][m]-coordinates[k][m])*(coordinates[j][m]-coordinates[k][m]);
								double d=Math.sqrt(dij)+Math.sqrt(djk)+Math.sqrt(djk);
								if(dmax<d) {
									dmax=d;
									imax=i;
									jmax=j;
									kmax=k;
								}
							}
						}
					}
				}
			}
		} */
		
		double b=0.0;
		int imax=0, jmax=0, kmax=0;
		for(int i=0; i<numAtoms; i++) {
			if(!atomLabel[i].trim().equals("H")) {	
				for(int j=1; j<numAtoms; j++) {
					if(!atomLabel[j].trim().equals("H")) {
						double dij=0.0;
						for(int m=0; m<3; m++) dij=dij+(coordinates[i][m]-coordinates[j][m])*(coordinates[i][m]-coordinates[j][m]);
						if(dij>b) {
							b=dij;
							imax=i;
							jmax=j;
						}
					}
				}
			}
		}
		b=Math.sqrt(b);
		double h=0.0;
		for(int k=0; k<numAtoms; k++) {
			if(!atomLabel[k].trim().equals("H")) {
				double a=0.0;
				double c=0.0;
				for(int m=0; m<3; m++) a=a+(coordinates[imax][m]-coordinates[k][m])*(coordinates[imax][m]-coordinates[k][m]);
				for(int m=0; m<3; m++) c=c+(coordinates[jmax][m]-coordinates[k][m])*(coordinates[jmax][m]-coordinates[k][m]);
				a=Math.sqrt(a);
				c=Math.sqrt(c);
				double hk=0.5*Math.sqrt((a+b+c)*(b+c-a)*(a-b+c)*(a+b-c))/b;
				if(hk>h) {
					h=hk;
					kmax=k;
				}
			}
		}
		
		
		try {
			PrintfWriter rotFile = new PrintfWriter(new FileWriter(fileName));
//			rotFile.println("xyzfiles");
//			rotFile.println(fp.trim()+na.trim()+".xyz");
//			rotFile.println(fp.trim()+nb.trim()+".xyz");
			rotFile.println("fragments");
			rotFile.print(fp.trim()+nb.trim());
			rotFile.println(" "+fp.trim()+na.trim());
			rotFile.println("target");
			rotFile.printf("%4d",imax+1); rotFile.printf("%4d",imax+1); rotFile.println();
			rotFile.printf("%4d",jmax+1); rotFile.printf("%4d",jmax+1); rotFile.println();
			rotFile.printf("%4d",kmax+1); rotFile.printf("%4d",kmax+1); rotFile.println();
			rotFile.println("states");
			for(int i=0; i<numStates; i++) {
				Integer index = ndxStateList[stateset][i];
				rotFile.print(stateNames[index].trim());	
                if(i!=numStates+1) rotFile.print(" ");
			}
			rotFile.println();
			rotFile.close();
		} catch(IOException ei) {
			return false;
		}
		return true;
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
						convergedName=projectRoot+fragmentNames[frag]+"_DFT.xyz";
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
	
	public void write_NWChem_DFT(String nameF, String nameP, Integer mult, Integer ranks, Double[] rt, Integer memory, String acc, String jobName, String limit) {
//		public void write_NWChem_DFT(Integer frag, Integer ranks, Integer stat) {
		String fileName = nameP.trim()+".xyz";
		String slurmName;
		String lsfName;
		fragmentName=nameP.trim();
		if(!read_XYZ()) {
			System.out.println("write_NWChem_DFT Reading "+fileName+" failed");
			System.exit(0);
				};
		for(int k=0; k<6; k++) RandT[k]=rt[k];
		rotate_AND_translate(RandT);
		fileName = nameP.trim()+"_DFT.nw";
		try {
			PrintfWriter nwFile = new PrintfWriter(new FileWriter(fileName));
			nwFile.println("start "+nameP);
			nwFile.println("echo");
			nwFile.println("basis \"ao basis\" print");
			nwFile.println("* library \"def2-tzvp\"");
			nwFile.println("end");
			nwFile.println("geometry units angstrom autosym");
		    for(int i=0; i<numAtoms; i++) {
		    	nwFile.print(atomLabel[i]+" "); 
				nwFile.printf("%16.8f",coordinates[i][0]);
				nwFile.printf("%16.8f",coordinates[i][1]);
				nwFile.printf("%16.8f",coordinates[i][2]);
				nwFile.println();
		    }
		    nwFile.println("end");
			nwFile.println("dft");
			nwFile.println(" xc b3lyp");
			if(mult>1) nwFile.println(" mult "+mult);
			nwFile.println("end");
			nwFile.println("driver");
			nwFile.println("end");
			nwFile.println("task dft optimize");
			nwFile.close();
		} catch(IOException e) {
		}
		fileName = nameP.trim()+"_DFT.run";
		slurmName = nameP.trim()+"_DFT.slurm";
		lsfName = nameP.trim()+"_DFT.lsf";
		try {
			PrintfWriter runFile = new PrintfWriter(new FileWriter(fileName));
			PrintfWriter slurmFile = new PrintfWriter(new FileWriter(slurmName));
			PrintfWriter lsfFile = new PrintfWriter(new FileWriter(lsfName));
			runFile.println("#!/usr/bin/tcsh");
			slurmFile.println("#!/usr/bin/tcsh");
			lsfFile.println("#!/usr/bin/tcsh");
			
			runFile.println("setenv MOLCAS_NPROCS "+ranks);
			runFile.println("setenv MOLCAS_MEM "+memory);

			slurmFile.println("#SBATCH --account="+acc.trim());
			slurmFile.println("#SBATCH --nodes=1");
			slurmFile.println("#SBATCH --tasks-per-node="+ranks);
			slurmFile.println("#SBATCH --job-name="+jobName.trim());
			slurmFile.println("#SBATCH --time="+limit.trim());
			slurmFile.println();
			slurmFile.println("module purge");
			slurmFile.println("module load intel openmpi");
			slurmFile.println();
			slurmFile.println("#source required files");
			slurmFile.println();
			slurmFile.println("#cd to job directory");
			slurmFile.println();
			
			lsfFile.println("#BSUB -P "+acc.trim());
			lsfFile.println("#BSUB -J "+jobName.trim());
			lsfFile.println("#BSUB -W "+limit.trim());
			lsfFile.println("#BSUB -q batch ");
			lsfFile.println("#BSUB -o job.o%J");
			lsfFile.println("#BSUB -e job.e%J");
			lsfFile.println("#BSUB -alloc_flags \"gpumps\"");
			lsfFile.println("#BSUB -nnodes "+ranks);
			lsfFile.println();
			lsfFile.println("#module purge");
			lsfFile.println("#module load intel openmpi");
			lsfFile.println();
			lsfFile.println("#source required files");
			lsfFile.println();
			lsfFile.println("#cd to job directory");
			lsfFile.println();
			
		    runFile.println("mpirun -n "+ranks+" nwchem "+nameP.trim()+"_DFT > "+nameP.trim()+"_DFT.nwout");
		    slurmFile.println("mpirun -n "+ranks+" nwchem "+nameP.trim()+"_DFT > "+nameP.trim()+"_DFT.nwout");
			lsfFile.println("jsrun -a1 -g1 -r6 -c7 -d cyclic ./nwchem "+nameP.trim()+"_DFT > "+nameP.trim()+"_DFT.nwout");
		    
			runFile.close();
			slurmFile.close();
			lsfFile.close();
		} catch(IOException e) {
		}	
	}
	
	public Boolean write_Molcas_Int(String nameF, String nameP, String bs, String ct, Integer ch, Double fx, Double fy, Double fz) {
		String fileName = nameP.trim()+".xyz";
		fragmentName=nameP.trim();
		String rootName=nameP.trim();
		if(!read_XYZ()) System.exit(0);
		fileName = nameP+"_INT.input";
		String previous;
		create_Basis(bs,ct);
		Integer ndx=0;
		try {
			PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
			inputFile.println("&seward");
			if(ch==0) inputFile.println("high cholesky");
			if(ch==1) inputFile.println("medium cholesky");
			if(ch==2) inputFile.println("low cholesky");
		    previous=" ";
		    Integer count=0;
		    for(int i=0; i<numAtoms; i++) {
		    	if(!atomLabel[i].trim().equals(previous)) {
		    		previous=atomLabel[i].trim();
		    		if(i>0) inputFile.println("end basis set");
		    		inputFile.println("basis set");
		    		ndx=basis_Index(previous);
		    		inputFile.println(basisSet[ndx]);
//		    		if(previous.equals("H")) inputFile.println("h.ano-s...3s2p.");
//		    		if(previous.equals("C")) inputFile.println("c.ano-s...4s3p2d.");
//		    		if(previous.equals("N")) inputFile.println("n.ano-s...4s3p2d.");
//		    		if(previous.equals("O")) inputFile.println("o.ano-s...4s3p2d.");
//		    		if(previous.equals("F")) inputFile.println("f.ano-s...4s3p2d.");
//		    		if(previous.equals("Cl")) inputFile.println("cl.ano-s...5s4p3d.");
//		    		if(previous.equals("Cu")) inputFile.println("cu.ano-s...7s5p4d3f.");
//		    		if(previous.equals("Br")) inputFile.println("br.ano-s...6s5p4d.");
		    		count=0;
		    	}
		    	count++;
		    	inputFile.print(atomLabel[i].trim()+count+" "); 
				inputFile.printf("%16.8f",coordinates[i][0]);
				inputFile.printf("%16.8f",coordinates[i][1]);
				inputFile.printf("%16.8f",coordinates[i][2]);
				inputFile.println(" Angstrom");
				previous=atomLabel[i].trim();
		    }
			inputFile.println("end basis set");
	    	inputFile.println();
		    if(fx!=0.0 || fy!=0.0 || fz!=0.0) {
		    	inputFile.println("&ffpt &end");
		    	inputFile.println("title");
		    	inputFile.println(" Add small electric field in x direction");
		    	inputFile.println("dipo");
		    	if(fx!=0) inputFile.println("x "+fx);
		    	if(fy!=0) inputFile.println("y "+fy);
		    	if(fz!=0) inputFile.println("z "+fz);
		    	inputFile.println("end of input");
		    	inputFile.println();
		    }
			inputFile.println(">>> COPY "+rootName.trim()+".OneInt $CurrDir/"+rootName.trim()+".oneint");
			inputFile.println(">>> COPY "+rootName.trim()+".RunFile $CurrDir/"+rootName.trim()+".runfil");
		    inputFile.close();
		    return true;
		} catch(IOException e) {
			return false;
		}
	}

	public Boolean write_Molcas_SCF(String nameF, String nameP, Integer mult, Integer chrg) {
		String fileName = nameP.trim()+".xyz";
		fragmentName=nameP.trim();
		String rootName=nameP.trim();
		if(!read_XYZ()) System.exit(0);
		fileName = nameP+"_SCF.input";
		try {
			PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
			inputFile.println("&scf");
			if(chrg!=0) {
				if(mult==1) { inputFile.println(" charge "); inputFile.println("  "+chrg);}
			} else {
				if(mult==2 || mult==4) {
					inputFile.println(" charge");inputFile.println("  1");
//					inputFile.println(" uhf");
//					inputFile.println(" spin "+mult);
				}
			}
			inputFile.close();
			return true;
		} catch(IOException e) {
			return false;
		}
	}
	
	public void Molcas_numAlt(String nameF,String nameP) {
		String fileName=nameP.trim()+".alter";
		String card;
		numAlt=0;
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			card=br.readLine();
			numAlt=Integer.valueOf(card.substring(0,6).trim());
			
			for(int i=0; i<numAlt; i++) {
				card=br.readLine();
				alter[i][0]=Integer.valueOf(card.substring(0,6).trim());
				alter[i][1]=Integer.valueOf(card.substring(6,12).trim());
			}
			br.close();
			return;
		} catch(IOException ef) {
			return;
		}
	}
	
	public Boolean write_Molcas_CASSCF(String nameF, String nameP, String nameS, Boolean withCASPT2, Integer numElec, Integer numCASe, Integer numCASo, Boolean withAlter, Double ipea) {
		Molcas_numAlt(nameF, nameP);
		String fileName=nameP.trim()+"_"+nameS.trim()+".input";
		String rootName=nameP.trim();
//		String fileName = projectRoot.trim()+fragmentNames[frag].trim()+"_"+stateNames[state].trim()+".input";
//		String rootName=projectRoot.trim()+fragmentNames[frag].trim();
		String ext = "_"+nameS.trim();
		altDone=!withAlter;
		try {
			Integer Inact = (numElec - numCASe)/2;
			if(nameS.trim().equals("S0")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				if(!altDone && numAlt>0) {
					altDone=true;
					inputFile.println("alter");
					inputFile.println(" "+numAlt);
					for(int i=0; i<numAlt; i++) {
						inputFile.println(" 1 "+alter[i][0]+" "+alter[i][1]);
					}
				}
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
//				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".orb");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				inputFile.println("&grid_it");
				inputFile.println("name=S0");
				inputFile.println("select");
				inputFile.println("1:"+(Inact+1)+"-"+(Inact+numCASo));
				inputFile.println("dense");
				inputFile.println(">>> COPY "+rootName.trim()+".S0.lus $CurrDir/"+rootName.trim()+ext.trim()+".lus");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("maxiter = 40");
					inputFile.println("ipea = "+ipea);
				}
				inputFile.close();
				return true;
			} else if(nameS.trim().equals("S1")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				if(!altDone && numAlt>0) {
					altDone=true;
					inputFile.println("alter");
					inputFile.println(" "+numAlt);
					for(int i=0; i<numAlt; i++) {
						inputFile.println(" 1 "+alter[i][0]+" "+alter[i][1]);
					}
				}
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
//				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.2 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.2 $CurrDir/"+rootName.trim()+ext.trim()+".orb");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.2 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				inputFile.println("&grid_it");
				inputFile.println("name=S1");
				inputFile.println("select");
				inputFile.println("1:"+(Inact+1)+"-"+(Inact+numCASo));
				inputFile.println("dense");
				inputFile.println(">>> COPY "+rootName.trim()+".S1.lus $CurrDir/"+rootName.trim()+ext.trim()+".lus");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("Multistate= 2 1 2");
					inputFile.println("maxiter = 30");
					inputFile.println("ipea = "+ipea);
				}
				inputFile.close();
				return true;
			} else if(nameS.trim().equals("S2")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				if(!altDone && numAlt>0) {
					altDone=true;
					inputFile.println("alter");
					inputFile.println(" "+numAlt);
					for(int i=0; i<numAlt; i++) {
						inputFile.println(" 1 "+alter[i][0]+" "+alter[i][1]);
					}
				}
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
//				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.3 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.3 $CurrDir/"+rootName.trim()+ext.trim()+".orb");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.3 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				inputFile.println("&grid_it");
				inputFile.println("name=S2");
				inputFile.println("select");
				inputFile.println("1:"+(Inact+1)+"-"+(Inact+numCASo));
				inputFile.println("dense");
				inputFile.println(">>> COPY "+rootName.trim()+".S2.lus $CurrDir/"+rootName.trim()+ext.trim()+".lus");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("Multistate= 2 1 2");
					inputFile.println("maxiter = 30");
					inputFile.println("ipea = "+ipea);
				}
				inputFile.close();
				return true;
			} else if(nameS.trim().equals("D0")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				if(!altDone && numAlt>0) {
					altDone=true;
					inputFile.println("alter");
					inputFile.println(" "+numAlt);
					for(int i=0; i<numAlt; i++) {
						inputFile.println(" 1 "+alter[i][0]+" "+alter[i][1]);
					}
				}
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
//				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".orb");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				inputFile.println("&grid_it");
				inputFile.println("name=D0");
				inputFile.println("select");
				inputFile.println("1:"+(Inact+1)+"-"+(Inact+numCASo));
				inputFile.println("dense");
				inputFile.println(">>> COPY "+rootName.trim()+".D0.lus $CurrDir/"+rootName.trim()+ext.trim()+".lus");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("maxiter = 30");
					inputFile.println("ipea = "+ipea);
				}
				inputFile.close();
				return true;
			} else if(nameS.trim().equals("D1")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				if(!altDone && numAlt>0) {
					altDone=true;
					inputFile.println("alter");
					inputFile.println(" "+numAlt);
					for(int i=0; i<numAlt; i++) {
						inputFile.println(" 1 "+alter[i][0]+" "+alter[i][1]);
					}
				}
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
//				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.2 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.2 $CurrDir/"+rootName.trim()+ext.trim()+".orb");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.2 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				inputFile.println("&grid_it");
				inputFile.println("name=D1");
				inputFile.println("select");
				inputFile.println("1:"+(Inact+1)+"-"+(Inact+numCASo));
				inputFile.println("dense");
				inputFile.println(">>> COPY "+rootName.trim()+".D1.lus $CurrDir/"+rootName.trim()+ext.trim()+".lus");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("Multistate= 2 1 2");
					inputFile.println("maxiter = 30");
					inputFile.println("ipea = "+ipea);
				}
				inputFile.close();
				return true;
			} else if(nameS.trim().equals("T1")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				if(!altDone && numAlt>0) {
					altDone=true;
					inputFile.println("alter");
					inputFile.println(" "+numAlt);
					for(int i=0; i<numAlt; i++) {
						inputFile.println(" 1 "+alter[i][0]+" "+alter[i][1]);
					}
				}
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
//				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".orb");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				inputFile.println("&grid_it");
				inputFile.println("name=T1");
				inputFile.println("select");
				inputFile.println("1:"+(Inact+1)+"-"+(Inact+numCASo));
				inputFile.println("dense");
				inputFile.println(">>> COPY "+rootName.trim()+".T1.lus $CurrDir/"+rootName.trim()+ext.trim()+".lus");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("maxiter = 30");
					inputFile.println("ipea = "+ipea);
				}
				inputFile.close();
				return true;
			} else if(nameS.trim().equals("T2")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				if(!altDone && numAlt>0) {
					altDone=true;
					inputFile.println("alter");
					inputFile.println(" "+numAlt);
					for(int i=0; i<numAlt; i++) {
						inputFile.println(" 1 "+alter[i][0]+" "+alter[i][1]);
					}
				}
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
//				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.2 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.2 $CurrDir/"+rootName.trim()+ext.trim()+".orb");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.2 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				inputFile.println("&grid_it");
				inputFile.println("name=T2");
				inputFile.println("select");
				inputFile.println("1:"+(Inact+1)+"-"+(Inact+numCASo));
				inputFile.println("dense");
				inputFile.println(">>> COPY "+rootName.trim()+".T2.lus $CurrDir/"+rootName.trim()+ext.trim()+".lus");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("Multistate= 2 1 2");
					inputFile.println("maxiter = 30");
					inputFile.println("ipea = "+ipea);
				}
				inputFile.close();
				return true;
			} else if(nameS.trim().equals("S+")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				if(!altDone && numAlt>0) {
					altDone=true;
					inputFile.println("alter");
					inputFile.println(" "+numAlt);
					for(int i=0; i<numAlt; i++) {
						inputFile.println(" 1 "+alter[i][0]+" "+alter[i][1]);
					}
				}
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
//				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".orb");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				inputFile.println("&grid_it");
				inputFile.println("name=S+");
				inputFile.println("select");
				inputFile.println("1:"+(Inact+1)+"-"+(Inact+numCASo));
				inputFile.println("dense");
				inputFile.println(">>> COPY "+rootName.trim()+".S+.lus $CurrDir/"+rootName.trim()+ext.trim()+".lus");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("maxiter = 30");
					inputFile.println("ipea = "+ipea);
				}
				inputFile.close();
				return true;
			} else if(nameS.trim().equals("D+")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				if(!altDone && numAlt>0) {
					altDone=true;
					inputFile.println("alter");
					inputFile.println(" "+numAlt);
					for(int i=0; i<numAlt; i++) {
						inputFile.println(" 1 "+alter[i][0]+" "+alter[i][1]);
					}
				}
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
//				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".orb");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				inputFile.println("&grid_it");
				inputFile.println("name=D+");
				inputFile.println("select");
				inputFile.println("1:"+(Inact+1)+"-"+(Inact+numCASo));
				inputFile.println("dense");
				inputFile.println(">>> COPY "+rootName.trim()+".D+.lus $CurrDir/"+rootName.trim()+ext.trim()+".lus");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("maxiter = 30");
					inputFile.println("ipea = "+ipea);
				}
				inputFile.close();
				return true;
			} else if(nameS.trim().equals("T+")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				if(!altDone && numAlt>0) {
					altDone=true;
					inputFile.println("alter");
					inputFile.println(" "+numAlt);
					for(int i=0; i<numAlt; i++) {
						inputFile.println(" 1 "+alter[i][0]+" "+alter[i][1]);
					}
				}
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
//				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".orb");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				inputFile.println("&grid_it");
				inputFile.println("name=T+");
				inputFile.println("select");
				inputFile.println("1:"+(Inact+1)+"-"+(Inact+numCASo));
				inputFile.println("dense");
				inputFile.println(">>> COPY "+rootName.trim()+".T+.lus $CurrDir/"+rootName.trim()+ext.trim()+".lus");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("maxiter = 30");
					inputFile.println("ipea = "+ipea);
				}
				inputFile.close();
				return true;
			} else if(nameS.trim().equals("S-")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				if(!altDone && numAlt>0) {
					altDone=true;
					inputFile.println("alter");
					inputFile.println(" "+numAlt);
					for(int i=0; i<numAlt; i++) {
						inputFile.println(" 1 "+alter[i][0]+" "+alter[i][1]);
					}
				}
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
//				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".orb");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				inputFile.println("&grid_it");
				inputFile.println("name=S-");
				inputFile.println("select");
				inputFile.println("1:"+(Inact+1)+"-"+(Inact+numCASo));
				inputFile.println("dense");
				inputFile.println(">>> COPY "+rootName.trim()+".S-.lus $CurrDir/"+rootName.trim()+ext.trim()+".lus");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("maxiter = 30");
					inputFile.println("ipea = "+ipea);
				}
				inputFile.close();
				return true;
			} else if(nameS.trim().equals("D-")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				if(!altDone && numAlt>0) {
					altDone=true;
					inputFile.println("alter");
					inputFile.println(" "+numAlt);
					for(int i=0; i<numAlt; i++) {
						inputFile.println(" 1 "+alter[i][0]+" "+alter[i][1]);
					}
				}
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
//				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".orb");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				inputFile.println("&grid_it");
				inputFile.println("name=D-");
				inputFile.println("select");
				inputFile.println("1:"+(Inact+1)+"-"+(Inact+numCASo));
				inputFile.println("dense");
				inputFile.println(">>> COPY "+rootName.trim()+".D-.lus $CurrDir/"+rootName.trim()+ext.trim()+".lus");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("maxiter = 30");
					inputFile.println("ipea = "+ipea);
				}
				inputFile.close();
				return true;
			} else if(nameS.trim().equals("T-")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				if(!altDone && numAlt>0) {
					altDone=true;
					inputFile.println("alter");
					inputFile.println(" "+numAlt);
					for(int i=0; i<numAlt; i++) {
						inputFile.println(" 1 "+alter[i][0]+" "+alter[i][1]);
					}
				}
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
//				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".orb");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				inputFile.println("&grid_it");
				inputFile.println("name=T-");
				inputFile.println("select");
				inputFile.println("1:"+(Inact+1)+"-"+(Inact+numCASo));
				inputFile.println("dense");
				inputFile.println(">>> COPY "+rootName.trim()+".T-.lus $CurrDir/"+rootName.trim()+ext.trim()+".lus");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("maxiter = 30");
					inputFile.println("ipea = "+ipea);
				}
				inputFile.close();
				return true;
			} else if(nameS.trim().equals("q1")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				if(!altDone && numAlt>0) {
					altDone=true;
					inputFile.println("alter");
					inputFile.println(" "+numAlt);
					for(int i=0; i<numAlt; i++) {
						inputFile.println(" 1 "+alter[i][0]+" "+alter[i][1]);
					}
				}
				inputFile.println("nactel");
				inputFile.println(" "+(numCASe));
				inputFile.println("spin");
				inputFile.println(" 4");
				inputFile.println("inactive");
				inputFile.println(" "+Inact);
				inputFile.println("ras2");
				inputFile.println(" "+numCASo);
				inputFile.println("prwf");
				inputFile.println("  0");
				inputFile.println("prsd");
//				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".orb");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				inputFile.println("&grid_it");
				inputFile.println("name=q1");
				inputFile.println("select");
				inputFile.println("1:"+(Inact+1)+"-"+(Inact+numCASo));
				inputFile.println("dense");
				inputFile.println(">>> COPY "+rootName.trim()+".q1.lus $CurrDir/"+rootName.trim()+ext.trim()+".lus");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("maxiter = 30");
					inputFile.println("ipea = "+ipea);
				}
				inputFile.close();
				return true;
			} else if(nameS.trim().equals("Q1")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println("&rasscf");
				if(!altDone && numAlt>0) {
					altDone=true;
					inputFile.println("alter");
					inputFile.println(" "+numAlt);
					for(int i=0; i<numAlt; i++) {
						inputFile.println(" 1 "+alter[i][0]+" "+alter[i][1]);
					}
				}
				inputFile.println("nactel");
				inputFile.println(" "+(numCASe));
				inputFile.println("spin");
				inputFile.println(" 5");
				inputFile.println("inactive");
				inputFile.println(" "+Inact);
				inputFile.println("ras2");
				inputFile.println(" "+numCASo);
				inputFile.println("prwf");
				inputFile.println("  0");
				inputFile.println("prsd");
//				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".orb");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				inputFile.println("&grid_it");
				inputFile.println("name=Q1");
				inputFile.println("select");
				inputFile.println("1:"+(Inact+1)+"-"+(Inact+numCASo));
				inputFile.println("dense");
				inputFile.println(">>> COPY "+rootName.trim()+".Q1.lus $CurrDir/"+rootName.trim()+ext.trim()+".lus");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("maxiter = 30");
					inputFile.println("ipea = "+ipea);
				}
				inputFile.close();
				return true;
			} else if(nameS.trim().equals("SQ1")) {
				PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
				inputFile.println(">>> COPY $CurrDir/"+rootName.trim()+"_Q1.orb INPORB");
				inputFile.println();
				inputFile.println("&rasscf");
				if(!altDone && numAlt>0) {
					altDone=true;
					inputFile.println("alter");
					inputFile.println(" "+numAlt);
					for(int i=0; i<numAlt; i++) {
						inputFile.println(" 1 "+alter[i][0]+" "+alter[i][1]);
					}
				}
				inputFile.println("nactel");
				inputFile.println(" "+(numCASe));
				inputFile.println("spin");
				inputFile.println(" 5");
				inputFile.println("inactive");
				inputFile.println(" "+Inact);
				inputFile.println("ras2");
				inputFile.println(" "+numCASo);
				inputFile.println("CIRoot");
				inputFile.println("  1 3");
				inputFile.println("  3");
				inputFile.println("lumorb");
				inputFile.println("cionly");
				inputFile.println("prwf");
				inputFile.println("  0");
				inputFile.println("prsd");
//				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".INPORB");
				inputFile.println(">>> COPY "+rootName.trim()+".RasOrb.1 $CurrDir/"+rootName.trim()+ext.trim()+".orb");
			    inputFile.println(">>> COPY "+rootName.trim()+".VecDet.1 $CurrDir/"+rootName.trim()+ext.trim()+".det");
				inputFile.println("&grid_it");
				inputFile.println("name=SQ1");
				inputFile.println("select");
				inputFile.println("1:"+(Inact+1)+"-"+(Inact+numCASo));
				inputFile.println("dense");
				inputFile.println(">>> COPY "+rootName.trim()+".SQ1.lus $CurrDir/"+rootName.trim()+ext.trim()+".lus");
				if(withCASPT2) {
					inputFile.println("&caspt2");
					inputFile.println("multistate = 3 1 2 3");
					inputFile.println("maxiter = 30");
					inputFile.println("ipea = "+ipea);
				}
				inputFile.close();
				return true;
			}
		} catch(IOException e) {
			return false;
		}
		return false;
	}
	
	public Boolean Molcas_SCF_Converged(Integer frag, Integer numCASe) {
		String fileName = fragmentName+fragmentNames[frag]+"_SCF.output";
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
				if(card.contains("++    Molecular itals:") && converged) {
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
		String fileName = fragmentName.trim()+fragmentNames[frag]+"_SCF.output";
		String card;
		Integer numOcc;
		Boolean converged = false;
		Double energy = 0.0;
		altDone=false;
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
				if(card.contains("++    Molecular oitals:") && converged) {
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
		String fileName = fragmentName.trim()+fragmentNames[frag]+"_"+stateNames[state]+".output";
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
					convergedCASSCF=true;
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
		String fileName = fragmentName.trim()+fragmentNames[frag]+"_"+stateNames[state]+".output";
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
		String fileName = fragmentName.trim()+fragmentNames[frag]+"_"+stateNames[state]+".output";
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

	public Boolean write_MEBF_XYZ(String p, String pn, Integer n, String[] fname, String[] frags, Double[][] randt) {
		String fileName = p+".xyz";
		Integer numberAtoms = 0;
		for(int j=0; j<n; j++) {
			for(int k=0; k<6; k++) RandT[k]=randt[j][k];
			initialize3(fname[j], frags[j], 2);
			numberAtoms=numberAtoms+numAtoms;
		}
		try {
			PrintfWriter xyzFile = new PrintfWriter(new FileWriter(fileName));
			xyzFile.printf("%6d",numberAtoms);
			xyzFile.println();
			xyzFile.println("Coordinates in Angstrom");
			for(int j=0; j<n; j++) {
				for(int k=0; k<6; k++) RandT[k]=randt[j][k];
				initialize3(fname[j], frags[j], 2);
				for(int i=0; i<numAtoms; i++) {
					xyzFile.print(atomLabel[i]); 
					xyzFile.printf("%16.8f",coordinates[i][0]);
					xyzFile.printf("%16.8f",coordinates[i][1]);
					xyzFile.printf("%16.8f",coordinates[i][2]);
					xyzFile.println();
				}
			}
			xyzFile.close();
			return true;
		} catch(IOException ei) {
			return false;
		}
	}
	
	public void write_Molcas_MEBF_One(String p, String pn, Integer n, String[] fname, String[] frags, Double[][] randt, String bs, String ct, Integer ch) {
		String fileName = p+"_ONE.input";
//		System.out.println("Constructing "+fileName);
		String previous;
		Integer atom=0;
		Integer maxElement = 45;
		Integer ndx=0;
		Integer[] counts = new Integer[maxElement];
		create_Basis(bs,ct);
		for(int i=0; i<maxElement; i++) counts[i]=0;
		try {
			PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
			inputFile.println("&seward");
			if(ch==0) inputFile.println("high cholesky");
			if(ch==1) inputFile.println("medium cholesky");
			if(ch==2) inputFile.println("low cholesky");
		    inputFile.println();
			
			for(int j=0; j<n; j++) {
				for(int k=0; k<6; k++) RandT[k]=randt[j][k];
				initialize3(fname[j], frags[j], 2);
			    previous=" ";
			    atom=0;
			    for(int i=0; i<numAtoms; i++) {
			    	if(!atomLabel[i].trim().equals(previous)) {
			    		previous=atomLabel[i].trim();
			    		if(i>0) inputFile.println("end basis set");
			    		inputFile.println("basis set");
			    		if(previous.equals("H")) atom=1;
			    		if(previous.equals("C")) atom=6;
			    		if(previous.equals("N")) atom=7;
			    		if(previous.equals("O")) atom=8;
			    		if(previous.equals("F")) atom=9;
			    		if(previous.equals("Cl")) atom=17;
			    		if(previous.equals("Cu")) atom=29;
			    		if(previous.equals("Br")) atom=35;
			    		for(int m=0; m<elements.length; m++) {
			    			if(atom==atomNumber[m]) inputFile.println(basisSet[m]);
			    		}
			    	}
			    	counts[atom]++;
			    	inputFile.print(atomLabel[i].trim()+counts[atom]+" "); 
					inputFile.printf("%16.8f",coordinates[i][0]);
					inputFile.printf("%16.8f",coordinates[i][1]);
					inputFile.printf("%16.8f",coordinates[i][2]);
					inputFile.println(" Angstrom");
					previous=atomLabel[i].trim();
			    }
			    inputFile.println("end basis set");
			    inputFile.println();
			}
			
		    inputFile.println("oneonly");
		    inputFile.println(">>> COPY "+p+".RunFile $CurrDir/"+p+".runfil");
		    inputFile.close();
		} catch(IOException e) {
		}
	}

	public void write_Molcas_MEBF_Two(String p, String pn, Integer n, String[] fname, String[] frags, Double[][] randt, String bs, String ct, Integer ch) {
		String fileName = p+"_TWO.input";
//		System.out.println("Constructing "+fileName);
		String previous;
		Integer atom=0;
		Integer maxElement = 45;
		Integer[] counts = new Integer[maxElement];
		create_Basis(bs,ct);
		for(int i=0; i<maxElement; i++) counts[i]=0;
		try {
			PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
			inputFile.println("&seward");
			if(ch==0) inputFile.println("high cholesky");
			if(ch==1) inputFile.println("medium cholesky");
			if(ch==2) inputFile.println("low cholesky");
		    inputFile.println();

			for(int j=0; j<n; j++) {
				for(int k=0; k<6; k++) RandT[k]=randt[j][k];
				initialize3(fname[j], frags[j], 2);
			    previous=" ";
			    atom=0;
			    for(int i=0; i<numAtoms; i++) {
			    	if(!atomLabel[i].trim().equals(previous)) {
			    		previous=atomLabel[i].trim();
			    		if(i>0) inputFile.println("end basis set");
			    		inputFile.println("basis set");
			    		if(previous.equals("H")) atom=1;
			    		if(previous.equals("C")) atom=6;
			    		if(previous.equals("N")) atom=7;
			    		if(previous.equals("O")) atom=8;
			    		if(previous.equals("F")) atom=9;
			    		if(previous.equals("Cl")) atom=17;
			    		if(previous.equals("Cu")) atom=29;
			    		if(previous.equals("Br")) atom=35;
			    		for(int m=0; m<elements.length; m++) {
			    			if(atom==atomNumber[m]) inputFile.println(basisSet[m]);
			    		}
			    	}
			    	counts[atom]++;
			    	inputFile.print(atomLabel[i].trim()+counts[atom]+" "); 
					inputFile.printf("%16.8f",coordinates[i][0]);
					inputFile.printf("%16.8f",coordinates[i][1]);
					inputFile.printf("%16.8f",coordinates[i][2]);
					inputFile.println(" Angstrom");
					previous=atomLabel[i].trim();
			    }
			    inputFile.println("end basis set");
			    inputFile.println();
			}

		    inputFile.println(">>> COPY $CurrDir/COMMONORB INPORB");
		    inputFile.println(">>> COPY $CurrDir/COMMONORB $CurrDir/"+p+".comorb");
		    inputFile.println();
		    inputFile.println("&motra");
		    inputFile.println("noorth");
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

	public void write_Molcas_MEBF_CB(String f, String p, Integer n, String[] frags, Integer[] fstat, Integer[] lenStateList, Integer[][] ndxStateList, Double thr) {
		String fileName = f+"_CB.input";
		try {
			PrintfWriter inputFile = new PrintfWriter(new FileWriter(fileName));
			inputFile.println("Project");
			inputFile.println(p);
			inputFile.println("Fragments");
			inputFile.printf("%3d",n); inputFile.println();
			for(int i=0; i<n; i++) inputFile.printf("%3d",lenStateList[fstat[i]]); inputFile.println();
			inputFile.println("Threshold");
			inputFile.println(" "+thr);
			inputFile.println("Labels");
			for(int i=0; i<n; i++) {
				inputFile.print(frags[i]);
				for(int j=0; j<lenStateList[fstat[i]]; j++) inputFile.print(" "+stateNames[ndxStateList[fstat[i]][j]].trim());
				inputFile.println();
			}
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
		g=Math.toRadians(rz);
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
		g=Math.toRadians(ry);
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
		g=Math.toRadians(rx);
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