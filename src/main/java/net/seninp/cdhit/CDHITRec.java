package net.seninp.cdhit;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

public class CDHITRec {

  // >Cluster 0
  // 0 22838nt, >ERR348569_k25_Locus_171_Trans... *
  // 1 967nt, >ERR348576_CL11246Contig1_1... at -/100.00%
  // 2 737nt, >ERR348576_CL14788Contig1_1... at +/99.86%
  // 3 2477nt, >ERR348576_CL353Contig1_1... at +/92.53%
  // 4 1734nt, >ERR348576_k25_Locus_32677_Tra... at -/100.00%
  // 5 346nt, >ERR348576_k25_Locus_48053_Tra... at +/100.00%
  // 6 1437nt, >ERR348576_k43_Locus_51503_Tra... at -/100.00%

  private long id;
  private CDHITClusterEntry centroid;
  private ArrayList<CDHITClusterEntry> entries;

  private CDHITRec() {
    super();
    entries = new ArrayList<CDHITClusterEntry>();
  }

  public static CDHITRec fromString(String buffer) {
    if (null == buffer || buffer.trim().isEmpty()) {
      return null;
    }
    // we parse the buffer line by line
    String[] lines = buffer.split(System.getProperty("line.separator"));
    CDHITRec res = new CDHITRec();
    res.setClusterId(toIdNumber(lines[0]));
    for (int i = 1; i < lines.length; i++) {
      String line = lines[i];
      if (line.endsWith("*")) {
        res.setCentroid(new CDHITClusterEntry(line));
      }
      else {
        res.addEntry(new CDHITClusterEntry(line));
      }
    }
    return res;
  }

  /**
   * Set the cluster's centroid.
   * 
   * @param cdhitClusterEntry the entry for centroid.
   */
  private void setCentroid(CDHITClusterEntry cdhitClusterEntry) {
    this.centroid = cdhitClusterEntry;
  }

  private void addEntry(CDHITClusterEntry e) {
    this.entries.add(e);
  }

  private void setClusterId(long id) {
    this.id = id;
  }

  private static long toIdNumber(String string) {
    return Long.valueOf(string.trim().split("\\s")[1]).longValue();
  }

  public long getId() {
    return this.id;
  }

  public CDHITClusterEntry getCentroid() {
    return this.centroid;
  }

  public int getSize() {
    return this.entries.size();
  }

  public List<CDHITClusterEntry> getEntries() {
    return this.entries;
  }

  public boolean containsLibs(String[] libs) {

    // make a list of libraries seen in da cluster
    HashSet<String> libraries = new HashSet<String>(this.entries.size());
    for (CDHITClusterEntry e : this.entries) {
      libraries.add(e.getLibraryName());
    }

    // check the centroid for membership
    if (!(libraries.contains(this.centroid.getLibraryName()))) {
      return false;
    }

    // check the membership
    for (String s : libs) {
      if (!(libraries.contains(s))) {
        return false;
      }
    }
    return true;

  }
}
