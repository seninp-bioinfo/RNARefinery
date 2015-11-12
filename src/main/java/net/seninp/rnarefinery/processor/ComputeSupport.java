package net.seninp.rnarefinery.processor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Map.Entry;
import java.util.Set;
import org.slf4j.LoggerFactory;
import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.Logger;
import net.seninp.rnarefinery.cdhit.CDHITClusterEntry;
import net.seninp.rnarefinery.cdhit.CDHITClusterProcessor;
import net.seninp.rnarefinery.cdhit.CDHITRec;

/**
 * Reads the clustering log, builds a CDHIT data structure, and prints output.
 * 
 * @author psenin
 * 
 */
public class ComputeSupport {

  // constants and formats
  //
  private static final String CR = "\n";
  private static final String COMMA = ",";
  private static final String SPACE = " ";
  private static final DecimalFormatSymbols otherSymbols = new DecimalFormatSymbols(Locale.US);
  private static DecimalFormat df = new DecimalFormat("0.00", otherSymbols);

  // data location
  //
  private static final String DATA_PREFIX = "/media/Stock/RNARefinery/preliminary0/";
  private static final String DATA_FNAME = "all_fpkm3_out.clstr";
  // private static final String DATA_PREFIX = "/media/Stock/RNARefinery/preliminary1/";
  // private static final String DATA_FNAME = "test2_all_fpkm3.cdhit_out.clstr";

  private static Logger consoleLogger;
  private static Level LOGGING_LEVEL = Level.INFO;

  static {
    consoleLogger = (Logger) LoggerFactory.getLogger(ComputeSupport.class);
    consoleLogger.setLevel(LOGGING_LEVEL);
  }

  /**
   * Main runnable.
   * 
   * @param args None are used.
   * @throws Exception If error occurs.
   */
  public static void main(String[] args) throws Exception {

    // read the data
    //
    consoleLogger.info("Reading " + DATA_FNAME + "...");
    HashMap<String, CDHITRec> clusters = CDHITClusterProcessor.readClusters(DATA_PREFIX
        + DATA_FNAME);
    consoleLogger.info("Read " + clusters.size() + " clusters in " + DATA_FNAME + " file.");

    // get the data stats
    //
    // clusters length dump
    BufferedWriter bw = new BufferedWriter(new FileWriter(new File(DATA_PREFIX
        + "clusters_length_samples.txt")));
    bw.write("cluster_id,run,length,support" + CR);
    for (Entry<String, CDHITRec> cluster : clusters.entrySet()) {
      String clusterName = cluster.getValue().getCentroid().getName();
      String runId = clusterName.substring(0, clusterName.indexOf("_"));
      bw.write(("\'" + cluster.getValue().getCentroid().getName() + "\'" + COMMA + "\'" + runId
          + "\'" + COMMA + cluster.getValue().getCentroid().getLength() + COMMA + getDistinctMasks(
            cluster.getValue()).size())
          + CR);
    }
    bw.close();
  }

  private static Set<String> getDistinctMasks(CDHITRec value) {
    HashSet<String> masks = new HashSet<String>(value.getEntries().size());
    for (CDHITClusterEntry entry : value.getEntries()) {
      String mask = entry.getName().substring(0, 9);
      masks.add(mask);
    }
    return masks;
  }
}
