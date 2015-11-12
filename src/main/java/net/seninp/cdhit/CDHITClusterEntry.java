package net.seninp.cdhit;

/**
 * Implements a single line of the cdhit cluster.
 * 
 * @author psenin
 * 
 */
public class CDHITClusterEntry {

  // 3 2477nt, >ERR348576_CL353Contig1_1... at +/92.53%
  private Long length;
  private String name;
  private String strand;
  private Double similarity;

  /**
   * Constructs an entry using that string.
   * 
   * @param line the input string.
   */
  public CDHITClusterEntry(String line) {

    String[] split = line.trim().split("[\\s\\,\\>]+");
    this.length = Long.valueOf(split[1].replace("nt", ""));
    this.name = split[2].replaceAll("\\.{3}", "");

    // treat the centroid properly
    // 0 22838nt, >ERR348569_k25_Locus_171_Trans... *\n
    if (line.trim().endsWith("*")) {
      this.strand = "+";
      this.similarity = 100.0D;
    }
    else {
      if (split[4].startsWith("+")) {
        this.strand = "+";
      }
      else if (split[4].startsWith("-")) {
        this.strand = "-";
      }
      this.similarity = Double.valueOf(split[4].replaceAll("[\\+\\-\\/\\%]", ""));
    }
  }

  /**
   * Get the entry length.
   * 
   * @return the entry length
   */
  public long getLength() {
    return length;
  }

  /**
   * Get the entry name.
   * 
   * @return the entry name.
   */
  public String getName() {
    return name;
  }

  /**
   * Get the entry strand.
   * 
   * @return the entry strand.
   */
  public String getStrand() {
    return strand;
  }

  /**
   * Get the similarity.
   * 
   * @return the similarity.
   */
  public double getSimilarity() {
    return similarity;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#hashCode()
   */
  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((length == null) ? 0 : length.hashCode());
    result = prime * result + ((name == null) ? 0 : name.hashCode());
    result = prime * result + ((similarity == null) ? 0 : similarity.hashCode());
    result = prime * result + ((strand == null) ? 0 : strand.hashCode());
    return result;
  }

  /*
   * (non-Javadoc)
   * 
   * @see java.lang.Object#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object obj) {
    if (this == obj) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    if (!(obj instanceof CDHITClusterEntry)) {
      return false;
    }
    CDHITClusterEntry other = (CDHITClusterEntry) obj;
    if (length == null) {
      if (other.length != null) {
        return false;
      }
    }
    else if (!length.equals(other.length)) {
      return false;
    }
    if (name == null) {
      if (other.name != null) {
        return false;
      }
    }
    else if (!name.equals(other.name)) {
      return false;
    }
    if (similarity == null) {
      if (other.similarity != null) {
        return false;
      }
    }
    else if (!similarity.equals(other.similarity)) {
      return false;
    }
    if (strand == null) {
      if (other.strand != null) {
        return false;
      }
    }
    else if (!strand.equals(other.strand)) {
      return false;
    }
    return true;
  }

  public String getLibraryName() {
    return this.name.split("_")[0];
  }

}
