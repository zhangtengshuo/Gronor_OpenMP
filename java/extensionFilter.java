package gronor;

import javax.swing.filechooser.FileFilter;
import java.io.File;

public class extensionFilter extends FileFilter
{

  private String extension;

  public extensionFilter(String ext)
  {
    extension=ext.toLowerCase();
  }
  public boolean accept(File file)
  {
    return (file.isDirectory() || file.getName().toLowerCase().endsWith(extension));
  }
  public String getDescription()
  {
    return " ";
  }
}
