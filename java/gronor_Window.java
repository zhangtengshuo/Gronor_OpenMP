package gronor;


import java.util.*;
import java.awt.*;
import java.io.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

public class gronor_Window extends JFrame implements ActionListener, ChangeListener, WindowListener, MouseListener{

	private static final long serialVersionUID = 4L;

	JMenuBar menubar = new JMenuBar();
	JMenuItem create, open, close, list, quit;
	
	gronor_Window(){

		super("GronOR Graphical User Interface");
	    super.setSize(800,200);
	    super.setLocation(50,50);
	    super.setBackground(Color.lightGray);
	    super.setForeground(Color.blue);
	    super.addWindowListener(this);
	    super.setFont(new Font("Serif",Font.BOLD,10));

	    this.setJMenuBar(menubar);

	    gronor_Project project;

	    JMenu select = new JMenu("File");
	    menubar.add(select);

	    select.add(open = new JMenuItem("Open Project"));
	    open.addActionListener(
	    			new ActionListener(){
	    				public void actionPerformed(ActionEvent e){
	    					JFileChooser chooser;
	    				    extensionFilter filter;
	    				    String prjFile;
	    				    String projectName;
	    				    JFrame dialogFrame;
	    					File prjF = new File("./");
	    				    FilenameFilter fnFilter = new FilenameFilter() {
	    				    	public boolean accept(File f, String name) {
	    				    		return name.endsWith(".prj");
	    				    	}
	    				    };
	    				    File[] files = prjF.listFiles(fnFilter);
	    				    if(files.length==1) {
	    				    	prjFile=files[0].getName();
	    				    } else {
	    				    	dialogFrame = new JFrame();
	    				    	dialogFrame.setSize(300,400);
	    				    	filter = new extensionFilter(".prj");
	    				    	chooser = new JFileChooser("./");
	    				    	chooser.setFileFilter(filter);
	    				    	chooser.showOpenDialog(dialogFrame);
	    				    	try {
	    				    	prjFile = chooser.getSelectedFile().toString();
	    				    	if(!prjFile.endsWith(".prj")) prjFile=prjFile.trim()+".prj";
	    				    	} catch(NullPointerException ee) {
	    				    		prjFile="Test.prj";
	    				    	}
	    				    };
	    					projectName=prjFile.substring(prjFile.lastIndexOf("/")+1,prjFile.indexOf(".prj"));
	    					gronor_Project project = new gronor_Project(projectName);
	    				} } );

	    select.add(close = new JMenuItem("Close Project"));
	    close.addActionListener(new ActionListener(){
	      public void actionPerformed(ActionEvent e){ 
		setVisible(false); System.exit(0); }});
	    
	    select.add(quit = new JMenuItem("Quit"));
	    quit.addActionListener(new ActionListener(){
	      public void actionPerformed(ActionEvent e){ 
		setVisible(false); System.exit(0); }});
	    
	    super.setVisible(true);
	    
		File prjF = new File("./");
	    FilenameFilter fnFilter = new FilenameFilter() {
	    	public boolean accept(File f, String name) {
	    		return name.endsWith(".prj");
	    	}
	    };
	    File[] files = prjF.listFiles(fnFilter);
	    if(files.length==1) {
		    String prjFile;
		    String projectName;
	    	prjFile=files[0].getName();
			projectName=prjFile.substring(prjFile.lastIndexOf("/")+1,prjFile.indexOf(".prj"));
			project = new gronor_Project(projectName);
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
