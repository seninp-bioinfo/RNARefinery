package net.seninp.cdhit;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;

public class CDHITClusterProcessor {

  private static final String CR = "\n";

  // private static final String COMMA = ", ";
  // private static final String SPACE = " ";

  /**
   * Disable constructor.
   */
  private CDHITClusterProcessor() {
    assert true;
  }

  /**
   * Formats an array to R sequence.
   * 
   * @param array
   * @return
   */
  public static String arrayToRSequence(Object[] array) {
    return "c(" + Arrays.toString(array).replaceAll("[\\[\\]\\s]+", "") + ")";
  }

  /**
   * Find the mean centroid length.
   * 
   * @param entry
   * @return
   */
  public static double meanCentroidLength(HashMap<String, CDHITRec> entry) {
    Double res = 0.;
    int counter = 0;
    for (Entry<String, CDHITRec> e : entry.entrySet()) {
      res = res + e.getValue().getCentroid().getLength();
      counter++;
    }
    return res / (double) counter;
  }

  /**
   * Get the longest centroid.
   * 
   * @param entry
   * @return
   */
  public static Integer maxCentroidLength(HashMap<String, CDHITRec> entry) {
    int res = Integer.MIN_VALUE;
    for (Entry<String, CDHITRec> e : entry.entrySet()) {
      long l = e.getValue().getCentroid().getLength();
      if (l > res) {
        res = (int) l;
      }
    }
    return res;
  }

  /**
   * Get the shortest centroid.
   * 
   * @param entry
   * @return
   */
  public static Integer minCentroidLength(HashMap<String, CDHITRec> entry) {
    int res = Integer.MAX_VALUE;
    for (Entry<String, CDHITRec> e : entry.entrySet()) {
      long l = e.getValue().getCentroid().getLength();
      if (l < res) {
        res = (int) l;
      }
    }
    return res;
  }

  /**
   * Reads the CDHIT clusters file.
   * 
   * @param fname the clusters file.
   * @return the clusters collection.
   * @throws Exception if error occurs.
   */
  public static HashMap<String, CDHITRec> readClusters(String fname) throws Exception {

    // the data structure: file -> clusters
    // cluster identified by a centroid
    // centroid is a string of form <length>@<sequence_id>
    HashMap<String, CDHITRec> clusters = new HashMap<String, CDHITRec>();

    //
    // down to business
    BufferedReader in = new BufferedReader(new FileReader(new File(fname)));

    // reading in
    String line = null;
    StringBuffer sb;
    line = in.readLine();
    while (null != line && line.startsWith(">")) {

      sb = new StringBuffer(2048);
      sb.append(line).append(CR);

      while (null != (line = in.readLine())) {
        if (line.startsWith(">")) {
          break;
        }
        sb.append(line).append(CR);
      }

      CDHITRec cluster = CDHITRec.fromString(sb.toString());

      String cname = cluster.getCentroid().getName() + "@" + cluster.getCentroid().getLength();

      if (clusters.containsKey(cname)) {
        in.close();
        throw new Exception("The cluster " + cname + " is already in the map!");
        // continue;
        // in.close();
        // throw new Exception("The cluster " + cname + " is already in the map!");
      }

      clusters.put(cname, cluster);

    }

    in.close();
    return clusters;
  }

}
