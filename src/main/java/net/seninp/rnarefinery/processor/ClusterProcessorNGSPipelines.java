package net.seninp.rnarefinery.processor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
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
public class ClusterProcessorNGSPipelines {

  // constants and formats
  //
  private static final String CR = "\n";
  private static final String COMMA = ",";
  // private static final String SPACE = " ";
  // private static final DecimalFormatSymbols otherSymbols = new DecimalFormatSymbols(Locale.US);
  // private static DecimalFormat df = new DecimalFormat("0.00", otherSymbols);

  // data location
  //
  private static final String DATA_PREFIX = "sandbox/";
  private static final String DATA_FNAME = "contigs_after_dna.cdhit.out.clstr";
  private static final Object QUOTE = "'";

  private static Logger consoleLogger;
  private static Level LOGGING_LEVEL = Level.INFO;

  static {
    consoleLogger = (Logger) LoggerFactory.getLogger(ClusterProcessorNGSPipelines.class);
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
    HashMap<String, CDHITRec> clusters = CDHITClusterProcessor
        .readClusters(DATA_PREFIX + DATA_FNAME);
    consoleLogger.info("Read " + clusters.size() + " clusters in " + DATA_FNAME + " file.");

    // get the data stats
    //
    // clusters length dump
    BufferedWriter bw = new BufferedWriter(
        new FileWriter(new File(DATA_PREFIX + "clusters_length_samples.txt")));
    for (Entry<String, CDHITRec> cluster : clusters.entrySet()) {
      bw.write(("\'" + cluster.getValue().getCentroid().getName() + "\'" + COMMA
          + cluster.getValue().getCentroid().getLength() + COMMA
          + getDistinctMasks(cluster.getValue()).size()) + CR);
    }
    bw.close();

    // sort the data by the cluster size
    //
    //
    ArrayList<CDHITRec> recs = new ArrayList<CDHITRec>();
    recs.addAll(clusters.values());
    Collections.sort(recs, new Comparator<CDHITRec>() {

      @Override
      public int compare(CDHITRec o1, CDHITRec o2) {
        int a = o1.getSize();
        int b = o2.getSize();
        if (a > b) {
          return 1;
        }
        else if (a < b) {
          return -1;
        }
        return 0;
      }

    });

    // extract clusters that embed ERR348568, ERR348580, SRR924585 only
    //
    bw = new BufferedWriter(new FileWriter(new File(DATA_PREFIX + "clusters_table.txt")));

    final String[] libs = { "ERR348568", "ERR348580", "SRR924585" };
    HashSet<String> libraries = new HashSet<String>(libs.length * 2);
    for (String s : libs) {
      libraries.add(s);
    }

    int limit = 117;
    int count = 0;
    for (int i = 0; i < recs.size(); i++) {

      // don't print more than needed
      if (count > limit) {
        break;
      }

      CDHITRec cluster = recs.get(i);

      // if all libs in da cluster -- print the table
      if (cluster.containsLibs(libs)
          && libraries.contains(cluster.getCentroid().getLibraryName())) {
        StringBuffer sb = new StringBuffer();
        String centroidField = "'" + cluster.getCentroid().getName() + "',";
        for (CDHITClusterEntry e : cluster.getEntries()) {
          if (libraries.contains(e.getLibraryName())
              && !(cluster.getCentroid().getName().equalsIgnoreCase(e.getName()))) {
            sb.append(centroidField);
            sb.append(QUOTE).append(e.getName()).append(QUOTE).append(COMMA);
            sb.append(QUOTE).append(e.getLibraryName()).append(QUOTE).append(COMMA);
            sb.append(e.getLength()).append(COMMA);
            sb.append(e.getSimilarity()).append(CR);
          }
        }
        bw.write(sb.toString());
        count++;
      }

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
