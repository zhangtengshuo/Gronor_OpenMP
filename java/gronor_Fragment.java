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
	
	String[] fragmentNames = new String[] {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M"};
    String[] stateNames = new String[] {"S0","S1","S2","D0","D1","T1","T2","S+","D+","T+","S-","D-","T-","q1","Q1","SQ1","S1+","S1-"};
    
    
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
	
	Integer numAlt=0;
	Integer[][] alter = new Integer[12][2];
	Integer numSup=0, numSupOrb=0;
	Integer[] supsym = new Integer[64];
	Boolean altDone=false;
	
	Integer numPz=0;
	Integer[] pzOrb = new Integer[100];

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
			runFile.println("# Copy the determinant file of reference fragment as it remains unchanged upon rotation");
			slurmFile.println("# Copy the determinant file of reference fragment as it remains unchanged upon rotation");
			lsfFile.println("# Copy the determinant file of reference fragment as it remains unchanged upon rotation");
			for(int i=0; i<numStates; i++) {
				Integer index = ndxStateList[stateset][i];
				runFile.println("cp "+fp.trim()+name2+"_"+stateNames[index]+".det "+fp.trim()+name1+"_"+stateNames[index]+".det");
				slurmFile.println("cp "+fp.trim()+name2+"_"+stateNames[index]+".det "+fp.trim()+name1+"_"+stateNames[index]+".det");
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
			
			if(nfrags>1) {
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
			}
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
			}
			runFile.println("cp "+p.trim()+".runfil RUNFIL");
			slurmFile.println("cp "+p.trim()+".runfil RUNFIL");
			lsfFile.println("cp "+p.trim()+".runfil RUNFIL");
			runFile.println("cp "+p.trim()+".oneint ONEINT");
			slurmFile.println("cp "+p.trim()+".oneint ONEINT");
			lsfFile.println("cp "+p.trim()+".oneint ONEINT");
			fullName = p.trim()+"_CB";
			runFile.println("gcommon < "+fullName.trim()+".input > "+fullName.trim()+".output");
			slurmFile.println("gcommon < "+fullName.trim()+".input > "+fullName.trim()+".output");
			lsfFile.println("gcommon < "+fullName.trim()+".input > "+fullName.trim()+".output");
			index=0;
			for(int i=0; i<nfrags; i++) {
				runFile.println("rm RUNFIL"+frags[i]);
				slurmFile.println("rm RUNFIL"+frags[i]);
				lsfFile.println("rm RUNFIL"+frags[i]);
			}
			
			fullName = p.trim()+"_CB";
			runFile.println("setenv DELETED ` grep \"Deleted orbitals in MOTRA\" "+fullName.trim()+".output | cut -b 30-34 `");
			slurmFile.println("setenv DELETED ` grep \"Deleted orbitals in MOTRA\" "+fullName.trim()+".output | cut -b 30-34 `");
			lsfFile.println("setenv DELETED ` grep \"Deleted orbitals in MOTRA\" "+fullName.trim()+".output | cut -b 30-34 `");
//			runFile.println("setenv FROZEN ` grep \"Frozen orbitals in MOTRA\" "+fullName.trim()+".output | cut -b 30-34 `");
//			slurmFile.println("setenv FROZEN ` grep \"Frozen orbitals in MOTRA\" "+fullName.trim()+".output | cut -b 30-34 `");
//			lsfFile.println("setenv FROZEN ` grep \"Frozen orbitals in MOTRA\" "+fullName.trim()+".output | cut -b 30-34 `");
			runFile.println("touch TRAINT");
			slurmFile.println("touch TRAINT");
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
			slurmFile.println("rdcho $MOLCAS_NPROCS > "+fullName.trim()+".output");
			slurmFile.println("rm _CHMOT1");
			lsfFile.println("rdcho $MOLCAS_NPROCS > "+fullName.trim()+".output");
			lsfFile.println("rm _CHMOT1");
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
		File f = new File(fileName);
		Boolean skip=false;
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
			} catch(IOException ei) {
				return false;
			};
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
			} catch(IOException ei) {
				return false;
			}
		}

		if(!write_POV()) return false;
		return true;

	}
	
	public Boolean write_POV() {
		
		double xMin = coordinates[0][0]-0.2;
		double yMin = coordinates[0][1]-0.2;
		double zMin = coordinates[0][2]-0.2;
		double xMax = coordinates[0][0]+0.2;
		double yMax = coordinates[0][1]+0.2;
		double zMax = coordinates[0][2]+0.2;
		for(int i=1; i<numAtoms; i++) {
			xMin=Math.min(coordinates[i][0],xMin);
			xMin=Math.min(coordinates[i][1],yMin);
			xMin=Math.min(coordinates[i][2],zMin);
			xMax=Math.max(coordinates[i][0],xMax);
			xMax=Math.max(coordinates[i][1],yMax);
			xMax=Math.max(coordinates[i][2],zMax);
		}
		
		String fileName = "colors.inc";
		File f = new File(fileName);
		if(!f.exists()) {
			try {
				PrintfWriter colorFile = new PrintfWriter(new FileWriter(fileName));
				colorFile.println("#declare Colors_Inc_Temp = version ;");
				colorFile.println("#version 2.0 ;");
				colorFile.println("#declare Gray05 = color red 0.05 green 0.05 blue 0.05 ;");
				colorFile.println("#declare Gray10 = color red 0.10 green 0.10 blue 0.10 ;");
				colorFile.println("#declare Gray15 = color red 0.15 green 0.15 blue 0.15 ;");
				colorFile.println("#declare Gray20 = color red 0.20 green 0.20 blue 0.20 ;");
				colorFile.println("#declare Gray25 = color red 0.25 green 0.25 blue 0.25 ;");
				colorFile.println("#declare Gray30 = color red 0.30 green 0.30 blue 0.30 ;");
				colorFile.println("#declare Gray35 = color red 0.35 green 0.35 blue 0.35 ;");
				colorFile.println("#declare Gray40 = color red 0.40 green 0.40 blue 0.40 ;");
				colorFile.println("#declare Gray45 = color red 0.45 green 0.45 blue 0.45 ;");
				colorFile.println("#declare Gray50 = color red 0.50 green 0.50 blue 0.50 ;");
				colorFile.println("#declare Gray55 = color red 0.55 green 0.55 blue 0.55 ;");
				colorFile.println("#declare Gray60 = color red 0.60 green 0.60 blue 0.60 ;");
				colorFile.println("#declare Gray65 = color red 0.65 green 0.65 blue 0.65 ;");
				colorFile.println("#declare Gray70 = color red 0.70 green 0.70 blue 0.70 ;");
				colorFile.println("#declare Gray75 = color red 0.75 green 0.75 blue 0.75 ;");
				colorFile.println("#declare Gray80 = color red 0.80 green 0.80 blue 0.80 ;");
				colorFile.println("#declare Gray85 = color red 0.85 green 0.85 blue 0.85 ;");
				colorFile.println("#declare Gray90 = color red 0.90 green 0.90 blue 0.90 ;");
				colorFile.println("#declare Gray95 = color red 0.95 green 0.95 blue 0.95 ;");
				colorFile.println("#declare Gray90 = color red 0.90 green 0.90 blue 0.90 ;");
				colorFile.println("#declare Gray95 = color red 0.95 green 0.95 blue 0.95 ;");
				colorFile.println("#declare DimGray = color red 0.329412 green 0.329412 blue 0.329412 ;");
				colorFile.println("#declare DimGrey = color red 0.329412 green 0.329412 blue 0.329412 ;");
				colorFile.println("#declare Gray = color red 0.752941 green 0.752941 blue 0.752941 ;");
				colorFile.println("#declare Grey = color red 0.752941 green 0.752941 blue 0.752941 ;");
				colorFile.println("#declare LightGray = color red 0.658824 green 0.658824 blue 0.658824 ;");
				colorFile.println("#declare LightGrey = color red 0.658824 green 0.658824 blue 0.658824 ;");
				colorFile.println("#declare VLightGrey = color red 0.80 green 0.80 blue 0.80 ;");
				colorFile.println("#declare White = color red 1.0 green 1.0 blue 1.0 ;");
				colorFile.println("#declare Red = color red 1.0 green 0.0 blue 0.0 ;");
				colorFile.println("#declare Green = color red 0.0 green 1.0 blue 0.0 ;");
				colorFile.println("#declare Blue = color red 0.0 green 0.0 blue 1.0 ;");
				colorFile.println("#declare Yellow = color red 1.0 green 1.0 blue 0.0 ;");
				colorFile.println("#declare Cyan = color red 0.0 green 1.0 blue 1.0 ;");
				colorFile.println("#declare Magenta = color red 1.0 green 0.0 blue 1.0 ;");
				colorFile.println("#declare Black = color red 0.0 green 0.0 blue 0.0 ;");
				colorFile.println("#declare Aquamarine = color red 0.439216 green 0.858824 blue 0.576471 ;");
				colorFile.println("#declare BlueViolet = color red 0.62352 green 0.372549 blue 0.623529 ;");
				colorFile.println("#declare Brown = color red 0.647059 green 0.164706 blue 0.164706 ;");
				colorFile.println("#declare CadetBlue = color red 0.372549 green 0.623529 blue 0.623529 ;");
				colorFile.println("#declare Coral = color red 1.0 green 0.498039 blue 0.0 ;");
				colorFile.println("#declare CornflowerBlue = color red 0.258824 green 0.258824 blue 0.435294 ;");
				colorFile.println("#declare DarkGreen = color red 0.184314 green 0.309804 blue 0.184314 ;");
				colorFile.println("#declare DarkOliveGreen = color red 0.309804 green 0.309804 blue 0.184314 ;");
				colorFile.println("#declare DarkOrchid = color red 0.6 green 0.196078 blue 0.8 ;");
				colorFile.println("#declare DarkSlateBlue = color red 0.419608 green 0.137255 blue 0.556863 ;");
				colorFile.println("#declare DarkSlateGray = color red 0.184314 green 0.309804 blue 0.309804 ;");
				colorFile.println("#declare DarkSlateGrey = color red 0.184314 green 0.309804 blue 0.309804 ;");
				colorFile.println("#declare DarkTurquoise = color red 0.439216 green 0.576471 blue 0.858824 ;");
				colorFile.println("#declare Firebrick = color red 0.556863 green 0.137255 blue 0.137255 ;");
				colorFile.println("#declare ForestGreen = color red 0.137255 green 0.556863 blue 0.137255 ;");
				colorFile.println("#declare Gold = color red 0.8 green 0.498039 blue 0.196078 ;");
				colorFile.println("#declare Goldenrod = color red 0.858824 green 0.858824 blue 0.439216 ;");
				colorFile.println("#declare GreenYellow = color red 0.576471 green 0.858824 blue 0.439216 ;");
				colorFile.println("#declare IndianRed = color red 0.309804 green 0.184314 blue 0.184314 ;");
				colorFile.println("#declare Khaki = color red 0.623529 green 0.623529 blue 0.372549 ;");
				colorFile.println("#declare LightBlue = color red 0.74902 green 0.847059 blue 0.847059 ;");
				colorFile.println("#declare LightSteelBlue = color red 0.560784 green 0.560784 blue 0.737255 ;");
				colorFile.println("#declare LimeGreen = color red 0.196078 green 0.8 blue 0.196078 ;");
				colorFile.println("#declare Maroon = color red 0.556863 green 0.137255 blue 0.419608 ;");
				colorFile.println("#declare MediumAquamarine = color red 0.196078 green 0.8 blue 0.6 ;");
				colorFile.println("#declare MediumBlue = color red 0.196078 green 0.196078 blue 0.8 ;");
				colorFile.println("#declare MediumForestGreen = color red 0.419608 green 0.556863 blue 0.137255 ;");
				colorFile.println("#declare MediumGoldenrod = color red 0.917647 green 0.917647 blue 0.678431 ;");
				colorFile.println("#declare MediumOrchid = color red 0.576471 green 0.439216 blue 0.858824 ;");
				colorFile.println("#declare MediumSeaGreen = color red 0.258824 green 0.435294 blue 0.258824 ;");
				colorFile.println("#declare MediumSlateBlue = color red 0.498039 green 1.0 blue 0.0 ;");
				colorFile.println("#declare MediumSpringGreen = color red 0.498039 green 1.0 blue 0.0 ;");
				colorFile.println("#declare MediumTurquoise = color red 0.439216 green 0.858824 blue 0.858824 ;");
				colorFile.println("#declare MediumVioletRed = color red 0.858824 green 0.439216 blue 0.576471 ;");
				colorFile.println("#declare MidnightBlue = color red 0.184314 green 0.184314 blue 0.309804 ;");
				colorFile.println("#declare Navy = color red 0.137255 green 0.137255 blue 0.556863 ;");
				colorFile.println("#declare NavyBlue = color red 0.137255 green 0.137255 blue 0.556863 ;");
				colorFile.println("#declare Orange = color red 1 green 0.5 blue 0.0 ;");
				colorFile.println("#declare OrangeRed = color red 1.0 green 0.498039 blue 0.0 ;");
				colorFile.println("#declare Orchid = color red 0.858824 green 0.439216 blue 0.858824 ;");
				colorFile.println("#declare PaleGreen = color red 0.560784 green 0.737255 blue 0.560784 ;");
				colorFile.println("#declare Pink = color red 0.737255 green 0.560784 blue 0.560784 ;");
				colorFile.println("#declare Plum = color red 0.917647 green 0.678431 blue 0.917647 ;");
				colorFile.println("#declare Salmon = color red 0.435294 green 0.258824 blue 0.258824 ;");
				colorFile.println("#declare SeaGreen = color red 0.137255 green 0.556863 blue 0.419608 ;");
				colorFile.println("#declare Sienna = color red 0.556863 green 0.419608 blue 0.137255 ;");
				colorFile.println("#declare SkyBlue = color red 0.196078 green 0.6 blue 0.8 ;");
				colorFile.println("#declare SlateBlue = color red 0.0 green 0.498039 blue 1.0 ;");
				colorFile.println("#declare SpringGreen = color red 0.0 green 1.0 blue 0.498039 ;");
				colorFile.println("#declare SteelBlue = color red 0.137255 green 0.419608 blue 0.556863 ;");
				colorFile.println("#declare Tan = color red 0.858824 green 0.576471 blue 0.439216 ;");
				colorFile.println("#declare Thistle = color red 0.847059 green 0.74902 blue 0.847059 ;");
				colorFile.println("#declare Turquoise = color red 0.678431 green 0.917647 blue 0.917647 ;");
				colorFile.println("#declare Violet = color red 0.309804 green 0.184314 blue 0.309804 ;");
				colorFile.println("#declare VioletRed = color red 0.8 green 0.196078 blue 0.6 ;");
				colorFile.println("#declare Wheat = color red 0.847059 green 0.847059 blue 0.74902 ;");
				colorFile.println("#declare YellowGreen = color red 0.6 green 0.8 blue 0.196078 ;");
				colorFile.println("#declare SummerSky = color red 0.22 green 0.69 blue 0.87 ;");
				colorFile.println("#declare RichBlue = color red 0.35 green 0.35 blue 0.67 ;");
				colorFile.println("#declare Brass = color red 0.71 green 0.65 blue 0.26 ;");
				colorFile.println("#declare Copper = color red 0.72 green 0.45 blue 0.20 ;");
				colorFile.println("#declare Bronze = color red 0.55 green 0.47 blue 0.14 ;");
				colorFile.println("#declare Bronze2 = color red 0.65 green 0.49 blue 0.24 ;");
				colorFile.println("#declare Silver = color red 0.90 green 0.91 blue 0.98 ;");
				colorFile.println("#declare BrightGold = color red 0.85 green 0.85 blue 0.10 ;");
				colorFile.println("#declare OldGold = color red 0.81 green 0.71 blue 0.23 ;");
				colorFile.println("#declare Feldspar = color red 0.82 green 0.57 blue 0.46 ;");
				colorFile.println("#declare Quartz = color red 0.85 green 0.85 blue 0.95 ;");
				colorFile.println("#declare Mica = color red 0.0 green 0.0 blue 0.0 ;");
				colorFile.println("#declare NeonPink = color red 1.00 green 0.43 blue 0.78 ;");
				colorFile.println("#declare DarkPurple = color red 0.53 green 0.12 blue 0.47 ;");
				colorFile.println("#declare NeonBlue = color red 0.30 green 0.30 blue 1.00 ;");
				colorFile.println("#declare CoolCopper = color red 0.85 green 0.53 blue 0.10 ;");
				colorFile.println("#declare MandarinOrange = color red 0.89 green 0.47 blue 0.20 ;");
				colorFile.println("#declare LightWood = color red 0.91 green 0.76 blue 0.65 ;");
				colorFile.println("#declare MediumWood = color red 0.65 green 0.50 blue 0.39 ;");
				colorFile.println("#declare DarkWood = color red 0.52 green 0.37 blue 0.26 ;");
				colorFile.println("#declare SpicyPink = color red 1.00 green 0.11 blue 0.68 ;");
				colorFile.println("#declare SemiSweetChoc = color red 0.42 green 0.26 blue 0.15 ;");
				colorFile.println("#declare BakersChoc = color red 0.36 green 0.20 blue 0.09 ;");
				colorFile.println("#declare Flesh = color red 0.96 green 0.80 blue 0.69 ;");
				colorFile.println("#declare NewTan = color red 0.92 green 0.78 blue 0.62 ;");
				colorFile.println("#declare NewMidnightBlue = color red 0.00 green 0.00 blue 0.61 ;");
				colorFile.println("#declare VeryDarkBrown = color red 0.35 green 0.16 blue 0.14 ;");
				colorFile.println("#declare DarkBrown = color red 0.36 green 0.25 blue 0.20 ;");
				colorFile.println("#declare DarkTan = color red 0.59 green 0.41 blue 0.31 ;");
				colorFile.println("#declare GreenCopper = color red 0.32 green 0.49 blue 0.46 ;");
				colorFile.println("#declare DkGreenCopper = color red 0.29 green 0.46 blue 0.43 ;");
				colorFile.println("#declare DustyRose = color red 0.52 green 0.39 blue 0.39 ;");
				colorFile.println("#declare HuntersGreen = color red 0.13 green 0.37 blue 0.31 ;");
				colorFile.println("#declare Scarlet = color red 0.55 green 0.09 blue 0.09 ;");
				colorFile.println("#declare Clear = color red 1.0 green 1.01 blue 1.0 filter 1.0 ;");
				colorFile.println("#declare Plane_Map = 0 ;");
				colorFile.println("#declare Sphere_Map = 1 ;");
				colorFile.println("#declare Cylinder_Map = 2 ;");
				colorFile.println("#declare Torus_Map = 5 ;");
				colorFile.println("#declare Bi   = 2 ;");
				colorFile.println("#declare Norm = 4 ;");
				colorFile.println("#version Colors_Inc_Temp ;");
				colorFile.close();
			} catch(IOException ei) {
				return false;
			}
		}

		fileName = "plane.inc";
		f = new File(fileName);
		if(!f.exists()) {
			try {
				PrintfWriter planeFile = new PrintfWriter(new FileWriter(fileName));
				planeFile.println("plane <0,1,0>, -1 pigment { White } }");
				planeFile.close();
			} catch(IOException ei) {
				return false;
			}
		}
				
		fileName = fragmentName+"-camera.inc";
		f = new File(fileName);
		try {
			PrintfWriter cameraFile = new PrintfWriter(new FileWriter(fileName));
			cameraFile.println("camera {");
			cameraFile.print("  location < 0.0 0.0 "); cameraFile.printf("%16.8f", -10.0*Math.max(Math.abs(xMax),Math.max(Math.abs(yMax),Math.abs(zMax)))); cameraFile.println();
			cameraFile.print("  look_at < "); cameraFile.printf("%16.8f",Math.abs(xMin)); cameraFile.print(" ");cameraFile.printf("%16.8f",Math.abs(yMin)); cameraFile.print(" ");cameraFile.printf("%16.8f",Math.abs(zMin)); cameraFile.println();
			cameraFile.println("  angle 20");
			cameraFile.println("}");
			cameraFile.println("light_source { <   0.0  0.0 -50.0 > color rgb < 1.0 1.0 1.0 > }");
			cameraFile.println("light_source { < -10.0 20.0 -10.0 > color rgb < 1.0 1.0 1.0 > }");
			cameraFile.println("light_source { <  10.0 20.0 -10.0 > color rgb < 1.0 1.0 1.0 > }");
			cameraFile.println("light_source { <   0.0 10.0  10.0 > color rgb < 1.0 1.0 1.0 > }");
			cameraFile.println("background { color rgb < 0.0 0.0 0.0 > }");
			cameraFile.close();
		} catch(IOException ei) {
			return false;
		}
		
		fileName = fragmentName+".pov";
		f = new File(fileName);
		Integer aNumber;
		Double aRadius;
		String aColor;
		Double cpk=1.5;
		Double scale=1.0;
		try {
			PrintfWriter povFile = new PrintfWriter(new FileWriter(fileName));
			povFile.println("#include "+fragmentName+"-camera.inc");
			povFile.println("#include colors.inc");
			povFile.println("#include plane.inc");
			for(int i=0; i<numAtoms; i++) {
				aNumber=getAtomNumber(atomLabel[i]);
				aRadius=getAtomRadius(aNumber);
				aColor=getAtomColor(aNumber);
				povFile.println("sphere {");
				povFile.print("  < "+scale*coordinates[i][0]+" "+scale*coordinates[i][1]+" "+scale*coordinates[i][2]+" > ");
				povFile.printf("%8.3f",scale*cpk*aRadius); povFile.println();
				povFile.println("  texture { pigment { "+aColor+" }");
				povFile.println("            finish { ambient 0.16");
				povFile.println("                     diffuse 0.48");
				povFile.println("                     phong 1.25");
				povFile.println("                     phong_size 200 }");
				povFile.println("}");
			}
			povFile.close();
		} catch(IOException ei) {
			return false;
		}
		
		return true;
	}
	
	public Integer getAtomNumber(String aL) {
        Integer number=0;
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
		if(aL.equals("Fr")) number=number+87;
		if(aL.equals("Ra")) number=number+88;
		if(aL.equals("Ac")) number=number+89;
		if(aL.equals("Th")) number=number+90;
		if(aL.equals("Pa")) number=number+91;
		if(aL.equals("U")) number=number+92;
		if(aL.equals("Np")) number=number+93;
		if(aL.equals("Pu")) number=number+94;
		if(aL.equals("Am")) number=number+95;
		if(aL.equals("Cm")) number=number+96;
		if(aL.equals("Bk")) number=number+97;
		if(aL.equals("Cf")) number=number+98;
		if(aL.equals("Es")) number=number+99;
		if(aL.equals("Fm")) number=number+100;
		if(aL.equals("Md")) number=number+101;
		if(aL.equals("No")) number=number+102;
		if(aL.equals("Lr")) number=number+103;
		if(aL.equals("Rf")) number=number+104;
		if(aL.equals("Db")) number=number+105;
		if(aL.equals("Sg")) number=number+106;
		if(aL.equals("Bh")) number=number+107;
		if(aL.equals("Hs")) number=number+108;
		if(aL.equals("Mt")) number=number+109;
		if(aL.equals("Ds")) number=number+110;
		if(aL.equals("Rg")) number=number+111;
		if(aL.equals("Cn")) number=number+112;
		if(aL.equals("Nh")) number=number+113;
		if(aL.equals("Fl")) number=number+114;
		if(aL.equals("Mc")) number=number+115;
		if(aL.equals("Lv")) number=number+116;
		if(aL.equals("Ts")) number=number+117;
		if(aL.equals("Og")) number=number+118;
		return number;
	}
	
	public Double getAtomRadius(Integer aN) {
		Double radius = 0.0;
		if(aN==   1) radius=0.35;
		if(aN==   2) radius=1.22;
		if(aN==   3) radius=1.23;
		if(aN==   4) radius=0.89;
		if(aN==   5) radius=0.88;
		if(aN==   6) radius=0.77;
		if(aN==   7) radius=0.70;
		if(aN==   8) radius=0.66;
		if(aN==   9) radius=0.58;
		if(aN==  10) radius=1.60;
		if(aN==  11) radius=1.40;
		if(aN==  12) radius=1.36;
		if(aN==  13) radius=1.25;
		if(aN==  14) radius=1.17;
		if(aN==  15) radius=1.10;
		if(aN==  16) radius=1.04;
		if(aN==  17) radius=0.99;
		if(aN==  18) radius=1.91;
		if(aN==  19) radius=2.03;
		if(aN==  20) radius=1.74;
		if(aN==  21) radius=1.44;
		if(aN==  22) radius=1.32;
		if(aN==  23) radius=1.22;
		if(aN==  24) radius=1.19;
		if(aN==  25) radius=1.17;
		if(aN==  26) radius=1.165;
		if(aN==  27) radius=1.16;
		if(aN==  28) radius=1.15;
		if(aN==  29) radius=1.17;
		if(aN==  30) radius=1.25;
		if(aN==  31) radius=1.25;
		if(aN==  32) radius=1.22;
		if(aN==  33) radius=1.21;
		if(aN==  34) radius=1.17;
		if(aN==  35) radius=1.14;
		if(aN==  36) radius=1.98;
		if(aN==  37) radius=2.22;
		if(aN==  38) radius=1.92;
		if(aN==  39) radius=1.62;
		if(aN==  40) radius=1.45;
		if(aN==  41) radius=1.34;
		if(aN==  42) radius=1.29;
		if(aN==  43) radius=1.27;
		if(aN==  44) radius=1.24;
		if(aN==  45) radius=1.25;
		if(aN==  46) radius=1.28;
		if(aN==  47) radius=1.34;
		if(aN==  48) radius=1.41;
		if(aN==  49) radius=1.50;
		if(aN==  50) radius=1.40;
		if(aN==  51) radius=1.41;
		if(aN==  52) radius=1.37;
		if(aN==  53) radius=1.33;
		if(aN==  54) radius=2.09;
		if(aN==  55) radius=2.35;
		if(aN==  56) radius=1.98;
		if(aN==  57) radius=1.69;
		if(aN==  58) radius=1.65;
		if(aN==  59) radius=1.65;
		if(aN==  60) radius=1.64;
		if(aN==  61) radius=1.65;
		if(aN==  62) radius=1.66;
		if(aN==  63) radius=1.65;
		if(aN==  64) radius=1.61;
		if(aN==  65) radius=1.59;
		if(aN==  66) radius=1.59;
		if(aN==  67) radius=1.58;
		if(aN==  68) radius=1.57;
		if(aN==  69) radius=1.56;
		if(aN==  70) radius=1.56;
		if(aN==  71) radius=1.56;
		if(aN==  72) radius=1.44;
		if(aN==  73) radius=1.34;
		if(aN==  74) radius=1.30;
		if(aN==  75) radius=1.28;
		if(aN==  76) radius=1.26;
		if(aN==  77) radius=1.26;
		if(aN==  78) radius=1.29;
		if(aN==  79) radius=1.34;
		if(aN==  80) radius=1.44;
		if(aN==  81) radius=1.55;
		if(aN==  82) radius=1.54;
		if(aN==  83) radius=1.52;
		if(aN==  84) radius=1.53;
		if(aN==  85) radius=1.50;
		if(aN==  86) radius=2.20;
		if(aN==  87) radius=3.24;
		if(aN==  88) radius=2.68;
		if(aN==  89) radius=2.25;
		if(aN==  90) radius=2.16;
		if(aN==  91) radius=1.93;
		if(aN==  92) radius=1.66;
		if(aN==  93) radius=1.57;
		if(aN==  94) radius=1.81;
		if(aN==  95) radius=2.21;
		if(aN==  96) radius=1.43;
		if(aN==  97) radius=1.42;
		if(aN==  98) radius=1.40;
		if(aN==  99) radius=1.39;
		if(aN== 100) radius=1.38;
		if(aN== 101) radius=1.37;
		if(aN== 102) radius=1.36;
		if(aN== 103) radius=1.34;
		if(aN== 104) radius=1.30;
		if(aN== 105) radius=1.30;
		return radius;
	}
	
	public String getAtomColor(Integer aN) {
		String aC="LightGray";
		if(aN==   1) aC="Blue";
		if(aN==   2) aC="Pink";
		if(aN==   3) aC="FireBrick";
		if(aN==   4) aC="SpicyPink";
		if(aN==   5) aC="Green";
		if(aN==   6) aC="LightGray";
		if(aN==   7) aC="SkyBlue";
		if(aN==   8) aC="Red";
		if(aN==   9) aC="GoldenRod";
		if(aN==  10) aC="SpicyPink";
		if(aN==  11) aC="Green";
		if(aN==  12) aC="DarkGreen";
		if(aN==  13) aC="Gray85";
		if(aN==  14) aC="GoldenRod";
		if(aN==  15) aC="Orange";
		if(aN==  16) aC="Yellow";
		if(aN==  17) aC="Green";
		return aC;
	}
	
	public Boolean write_rotharm_input(String fp, String na, String nb, Integer stateset, Integer[] lenStateList, Integer[][] ndxStateList) {
		String fileName = fp.trim()+na.trim()+"_rotate.input";
		Integer numStates = lenStateList[stateset];
		
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
		String fileName = nameP.trim()+".xyz";
		String slurmName;
		String lsfName;
		fragmentName=nameP.trim();
		if(!read_XYZ()) {
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
			inputFile.println(" pror");
			inputFile.println(" 1 2");
			if(chrg!=0) {
				if(mult==1) { inputFile.println(" charge "); inputFile.println("  "+chrg);}
			} else {
				if(mult==2 || mult==4) {
					inputFile.println(" charge");inputFile.println("  1");
				}
			}
			inputFile.close();
			return true;
		} catch(IOException e) {
			return false;
		}
	}
	
	public Boolean write_Molcas_CASSCF(String nameF, String nameP, String nameS, Boolean withCASPT2, Integer numElec, Integer numCASe, Integer numCASo, Boolean withAlter, Boolean withSup, Double ipea) {

		String fileName=nameP.trim()+"_"+nameS.trim()+".input";
		String rootName=nameP.trim();
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
				if(withSup && numSup>0) {
					inputFile.println("supsym");
					inputFile.println(" 1");
					inputFile.print(" "+numSup);
					for(int i=0; i<numSup; i++) {
						inputFile.print(" "+supsym[i]);
						}
					inputFile.println();
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
				if(withSup && numSup>0) {
					inputFile.println("supsym");
					inputFile.println(" 1");
					inputFile.print(" "+numSup);
					for(int i=0; i<numSup; i++) {
						inputFile.print(" "+supsym[i]);
						}
					inputFile.println();
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
				if(withSup && numSup>0) {
					inputFile.println("supsym");
					inputFile.println(" 1");
					inputFile.print(" "+numSup);
					for(int i=0; i<numSup; i++) {
						inputFile.print(" "+supsym[i]);
						}
					inputFile.println();
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
				if(withSup && numSup>0) {
					inputFile.println("supsym");
					inputFile.println(" 1");
					inputFile.print(" "+numSup);
					for(int i=0; i<numSup; i++) {
						inputFile.print(" "+supsym[i]);
						}
					inputFile.println();
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
				if(withSup && numSup>0) {
					inputFile.println("supsym");
					inputFile.println(" 1");
					inputFile.print(" "+numSup);
					for(int i=0; i<numSup; i++) {
						inputFile.print(" "+supsym[i]);
						}
					inputFile.println();
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
				if(withSup && numSup>0) {
					inputFile.println("supsym");
					inputFile.println(" 1");
					inputFile.print(" "+numSup);
					for(int i=0; i<numSup; i++) {
						inputFile.print(" "+supsym[i]);
						}
					inputFile.println();
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
				if(withSup && numSup>0) {
					inputFile.println("supsym");
					inputFile.println(" 1");
					inputFile.print(" "+numSup);
					for(int i=0; i<numSup; i++) {
						inputFile.print(" "+supsym[i]);
						}
					inputFile.println();
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
				if(withSup && numSup>0) {
					inputFile.println("supsym");
					inputFile.println(" 1");
					inputFile.print(" "+numSup);
					for(int i=0; i<numSup; i++) {
						inputFile.print(" "+supsym[i]);
						}
					inputFile.println();
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
				if(withSup && numSup>0) {
					inputFile.println("supsym");
					inputFile.println(" 1");
					inputFile.print(" "+numSup);
					for(int i=0; i<numSup; i++) {
						inputFile.print(" "+supsym[i]);
						}
					inputFile.println();
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
				if(withSup && numSup>0) {
					inputFile.println("supsym");
					inputFile.println(" 1");
					inputFile.print(" "+numSup);
					for(int i=0; i<numSup; i++) {
						inputFile.print(" "+supsym[i]);
						}
					inputFile.println();
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
				if(withSup && numSup>0) {
					inputFile.println("supsym");
					inputFile.println(" 1");
					inputFile.print(" "+numSup);
					for(int i=0; i<numSup; i++) {
						inputFile.print(" "+supsym[i]);
						}
					inputFile.println();
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
				if(withSup && numSup>0) {
					inputFile.println("supsym");
					inputFile.println(" 1");
					inputFile.print(" "+numSup);
					for(int i=0; i<numSup; i++) {
						inputFile.print(" "+supsym[i]);
						}
					inputFile.println();
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
				if(withSup && numSup>0) {
					inputFile.println("supsym");
					inputFile.println(" 1");
					inputFile.print(" "+numSup);
					for(int i=0; i<numSup; i++) {
						inputFile.print(" "+supsym[i]);
						}
					inputFile.println();
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
				if(withSup && numSup>0) {
					inputFile.println("supsym");
					inputFile.println(" 1");
					inputFile.print(" "+numSup);
					for(int i=0; i<numSup; i++) {
						inputFile.print(" "+supsym[i]);
						}
					inputFile.println();
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
				if(withSup && numSup>0) {
					inputFile.println("supsym");
					inputFile.println(" 1");
					inputFile.print(" "+numSup);
					for(int i=0; i<numSup; i++) {
						inputFile.print(" "+supsym[i]);
						}
					inputFile.println();
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
				if(withSup && numSup>0) {
					inputFile.println("supsym");
					inputFile.println(" 1");
					inputFile.print(" "+numSup);
					for(int i=0; i<numSup; i++) {
						inputFile.print(" "+supsym[i]);
						}
					inputFile.println();
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
			} else if(nameS.trim().equals("S1+")) {
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
				if(withSup && numSup>0) {
					inputFile.println("supsym");
					inputFile.println(" 1");
					inputFile.print(" "+numSup);
					for(int i=0; i<numSup; i++) {
						inputFile.print(" "+supsym[i]);
						}
					inputFile.println();
				}
				inputFile.println("nactel");
				inputFile.println(" "+(numCASe-1));
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
				
			} else if(nameS.trim().equals("S1-")) {
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
				if(withSup && numSup>0) {
					inputFile.println("supsym");
					inputFile.println(" 1");
					inputFile.print(" "+numSup);
					for(int i=0; i<numSup; i++) {
						inputFile.print(" "+supsym[i]);
						}
					inputFile.println();
				}
				inputFile.println("nactel");
				inputFile.println(" "+(numCASe+1));
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
				
			}
		} catch(IOException e) {
			return false;
		}
		return false;
	}
	
	public Boolean Molcas_SCF_Converged(String nameP, Integer frag, Integer numCASe, Integer numCASo) {
		String fileName = nameP+fragmentNames[frag]+"_SCF.output";
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
				if(card.contains("      Total number of orbitals")) {
					numOrb=Integer.valueOf(card.substring(33,38).trim());
				}
				if(card.contains("::    Total SCF energy")) {
					converged=true;
					energy=Double.valueOf(card.substring(51,68).trim()).doubleValue();
				}
				if(card.contains("++    Molecular orbitals:") && converged) {
					numBas=numOcc+numSec;
					Double[] occ = new Double[numOrb];
					Integer[] typ = new Integer[numOrb];
					Double[] coef = new Double[10];
					Double[] cmax = new Double[10];
					Integer[] imax = new Integer[10];
					Integer ilast;
					Integer count=0;
					Integer ndx=0;
					
					Boolean[] isPz = new Boolean[numOrb];
					
					numPz=0;
					while(!card.contains("Index")) card=br.readLine();
					while(!card.contains("--")) {
						card=br.readLine();
						if(card.length()>=5) {
							if(card.substring(0,5).trim().length()>0) {
								ndx=Integer.valueOf(card.substring(0,5).trim());
								if(ndx>0) isPz[(ndx-1)]=card.substring(45,47).trim().contains("pz");
							}
						}
					}
					
					numOrb=ndx;
					
					Integer n=numCASe/2, numpz=0;
					if(numCASe%2!=0) n++;
					
					Integer numiop=0;
					for(int i=0; i<numOcc-n; i++) {
						if(isPz[i]) numiop++;
					}
					Integer numaonp=0;
					for(int i=numOcc-n; i<numOcc; i++) {
						if(!isPz[i]) numaonp++;
					}
					Integer numaunp=0;
					for(int i=numOcc; i<numOcc+numCASo-n; i++) {
						if(!isPz[i]) numaunp++;
					}
					Integer numiup=0;
					for(int i=numOcc+numCASo-n; i<numOrb; i++) {
						if(isPz[i]) numiup++;
					}
					
					numAlt=0;
					for(int i=numOcc-n; i<numOcc; i++) {
						if(!isPz[i]) {
							for(int j=numOcc-n-1; j>=0; j--) {
								if(!isPz[i] && isPz[j]) {
									alter[numAlt][0]=(j+1);
									alter[numAlt][1]=(i+1);
									numAlt++;
									isPz[i]=true;
									isPz[j]=false;
									numiop--;
									numaonp--;
								}
							}
						}						
					}

					for(int i=numOcc; i<numOcc+n; i++) {
						if(!isPz[i]) {
							for(int j=numOcc+n; j<numOrb; j++) {
								if(!isPz[i] && isPz[j]) {
									alter[numAlt][0]=(i+1);
									alter[numAlt][1]=(j+1);
									numAlt++;
									isPz[i]=true;
									isPz[j]=false;
									numiup--;
									numaunp--;
								}
							}
						}
					}
					
					numSup=0;
					for(int i=0; i<numOcc-n; i++) {
						if(isPz[i]) {
							supsym[numSup]=i+1;
							numSup++;
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

	public Double Molcas_SCF(Integer frag, Integer numCASe, Integer numCASo) {
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
				if(card.contains("      Total number of orbitals")) {
					numOrb=Integer.valueOf(card.substring(33,38).trim());
				}
				if(card.contains("::    Total SCF energy")) {
					converged=true;
					energy=Double.valueOf(card.substring(51,68).trim()).doubleValue();
				}

				if(card.contains("++    Molecular orbitals:") && converged) {
					numBas=numOcc+numSec;
					Double[] occ = new Double[numOrb];
					Integer[] typ = new Integer[numOrb];
					Double[] coef = new Double[10];
					Double[] cmax = new Double[10];
					Integer[] imax = new Integer[10];
					Integer ilast;
					Integer count=0;
					Integer ndx=0;
					
					Boolean[] isPz = new Boolean[numOrb];
					
					numPz=0;
					while(!card.contains("Index")) card=br.readLine();
					while(!card.contains("--")) {
						card=br.readLine();
						if(card.length()>=5) {
							if(card.substring(0,5).trim().length()>0) {
								ndx=Integer.valueOf(card.substring(0,5).trim());
								if(ndx>0) isPz[(ndx-1)]=card.substring(45,47).trim().contains("pz");
							}
						}
					}
					
					numOrb=ndx;
					
					Integer n=numCASe/2, numpz=0;
					if(numCASe%2!=0) n++;
					
					Integer numiop=0;
					for(int i=0; i<numOcc-n; i++) {
						if(isPz[i]) numiop++;
					}
					Integer numaonp=0;
					for(int i=numOcc-n; i<numOcc; i++) {
						if(!isPz[i]) numaonp++;
					}
					Integer numaunp=0;
					for(int i=numOcc; i<numOcc+numCASo-n; i++) {
						if(!isPz[i]) numaunp++;
					}
					Integer numiup=0;
					for(int i=numOcc+numCASo-n; i<numOrb; i++) {
						if(isPz[i]) numiup++;
					}
					
					numAlt=0;
					for(int i=numOcc-n; i<numOcc; i++) {
						if(!isPz[i]) {
							for(int j=numOcc-n-1; j>=0; j--) {
								if(!isPz[i] && isPz[j]) {
									alter[numAlt][0]=(j+1);
									alter[numAlt][1]=(i+1);
									numAlt++;
									isPz[i]=true;
									isPz[j]=false;
									numiop--;
									numaonp--;
								}
							}
						}						
					}

					for(int i=numOcc+numCASo-1; i>=numOcc; i--) {
						if(!isPz[i]) {
							for(int j=numOcc+numCASo-n; j<numOrb; j++) {
								if(!isPz[i] && isPz[j]) {
									alter[numAlt][0]=(i+1);
									alter[numAlt][1]=(j+1);
									numAlt++;
									isPz[i]=true;
									isPz[j]=false;
									numiup--;
									numaunp--;
								}
							}
						}
					}
					
					numSup=0;
					for(int i=0; i<numOcc-n; i++) {
						if(isPz[i]) {
							supsym[numSup]=i+1;
							numSup++;
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
		} catch(IOException ei) {
			return false;
		}
		if(!write_MEBF_POV(p,n,fname,frags, randt)) return false;
		return true;
	}
	
	public Boolean write_MEBF_POV(String p, Integer n, String[] fname, String[] frags, Double[][] randt) {
		
		Integer numberAtoms = 0;
		for(int j=0; j<n; j++) {
			for(int k=0; k<6; k++) RandT[k]=randt[j][k];
			initialize3(fname[j], frags[j], 2);
			numberAtoms=numberAtoms+numAtoms;
		}
		
		Double xMin = coordinates[0][0]-0.2;
		Double yMin = coordinates[0][1]-0.2;
		Double zMin = coordinates[0][2]-0.2;
		Double xMax = coordinates[0][0]+0.2;
		Double yMax = coordinates[0][1]+0.2;
		Double zMax = coordinates[0][2]+0.2;
		
		for(int j=0; j<n; j++) {
			initialize3(fname[j], frags[j], 2);
			for(int i=0; i<numAtoms; i++) {
				xMin=Math.min(coordinates[i][0],xMin);
				xMin=Math.min(coordinates[i][1],yMin);
				xMin=Math.min(coordinates[i][2],zMin);
				xMax=Math.max(coordinates[i][0],xMax);
				xMax=Math.max(coordinates[i][1],yMax);
				xMax=Math.max(coordinates[i][2],zMax);
			}
		}
				
		String fileName = p+"-camera.inc";
		File f = new File(fileName);
		try {
			PrintfWriter cameraFile = new PrintfWriter(new FileWriter(fileName));
			cameraFile.println("camera {");
			cameraFile.print("  location < 0.0 0.0 "); cameraFile.printf("%16.8f", -10.0*Math.max(Math.abs(xMax),Math.max(Math.abs(yMax),Math.abs(zMax)))); cameraFile.println();
			cameraFile.print("  look_at < "); cameraFile.printf("%16.8f",Math.abs(xMin)); cameraFile.print(" ");cameraFile.printf("%16.8f",Math.abs(yMin)); cameraFile.print(" ");cameraFile.printf("%16.8f",Math.abs(zMin)); cameraFile.println();
			cameraFile.println("  angle 20");
			cameraFile.println("}");
			cameraFile.println("light_source { <   0.0  0.0 -50.0 > color rgb < 1.0 1.0 1.0 > }");
			cameraFile.println("light_source { < -10.0 20.0 -10.0 > color rgb < 1.0 1.0 1.0 > }");
			cameraFile.println("light_source { <  10.0 20.0 -10.0 > color rgb < 1.0 1.0 1.0 > }");
			cameraFile.println("light_source { <   0.0 10.0  10.0 > color rgb < 1.0 1.0 1.0 > }");
			cameraFile.println("background { color rgb < 0.0 0.0 0.0 > }");
			cameraFile.close();
		} catch(IOException ei) {
			return false;
		}
		
		fileName = p+".pov";
		f = new File(fileName);
		Integer aNumber;
		Double aRadius;
		String aColor;
		Double cpk=1.5;
		Double scale=1.0;
		try {
			PrintfWriter povFile = new PrintfWriter(new FileWriter(fileName));
			povFile.println("#include "+p+"-camera.inc");
			povFile.println("#include colors.inc");
			povFile.println("#include plane.inc");

			for(int j=0; j<n; j++) {
				initialize3(fname[j], frags[j], 2);
				for(int i=0; i<numAtoms; i++) {
					aNumber=getAtomNumber(atomLabel[i]);
					aRadius=getAtomRadius(aNumber);
					aColor=getAtomColor(aNumber);
					povFile.println("sphere {");
					povFile.print("  < "+scale*coordinates[i][0]+" "+scale*coordinates[i][1]+" "+scale*coordinates[i][2]+" > ");
					povFile.printf("%8.3f",scale*cpk*aRadius); povFile.println();
					povFile.println("  texture { pigment { "+aColor+" }");
					povFile.println("            finish { ambient 0.16");
					povFile.println("                     diffuse 0.48");
					povFile.println("                     phong 1.25");
					povFile.println("                     phong_size 200 }");
					povFile.println("}");
				}
			}
			povFile.close();
		} catch(IOException ei) {
			return false;
		}
		
		
		return true;
	}
	
	public void write_Molcas_MEBF_One(String p, String pn, Integer n, String[] fname, String[] frags, Double[][] randt, String bs, String ct, Integer ch) {
		String fileName = p+"_ONE.input";
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
		    inputFile.println(">>> COPY "+p+".OneInt $CurrDir/"+p+".oneint");
		    inputFile.close();
		} catch(IOException e) {
		}
	}

	public void write_Molcas_MEBF_Two(String p, String pn, Integer n, String[] fname, String[] frags, Double[][] randt, String bs, String ct, Integer ch, Integer nfr) {
		String fileName = p+"_TWO.input";
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
		    inputFile.println(" "+nfr);
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
		} catch(IOException e) {
		}
	}

	public void write_Molcas_MEBF_CB(String f, String p, Integer n, String[] frags, Integer[] fstat, Integer[] nfrz, Integer[] lenStateList, Integer[][] ndxStateList, Double thr) {
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
			inputFile.println("Frozen");
			for(int i=0; i<n; i++) {
				inputFile.printf("%3d",nfrz[i]);
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