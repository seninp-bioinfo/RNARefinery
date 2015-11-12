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
 * 
 */

// DROP TABLE IF EXISTS `clusters`;
// CREATE TABLE `clusters` (
// `centroid` varchar(256) NOT NULL,
// `member` varchar(256) NOT NULL,
// `similarity` float DEFAULT NULL,
// FULLTEXT KEY `centr_name` (`centroid`),
// FULLTEXT KEY `seq_name` (`member`)
// ) ENGINE=MyISAM DEFAULT CHARSET=latin1;
// LOCK TABLES `clusters` WRITE;
// ALTER TABLE `clusters` DISABLE KEYS;
// LOAD DATA LOCAL INFILE '/home/psenin/git/RNARefinery.git/sandbox/clusters_table.txt' INTO TABLE
// `clusters` COLUMNS TERMINATED BY '\t';
// ALTER TABLE `clusters` ENABLE KEYS;
// UNLOCK TABLES;
// alter table clusters add `centroid_run_id` int(11) after centroid;
// alter table clusters add `member_run_id` int(11) after `member`;
// update clusters set centroid_run_id=(SELECT sr.id from sra_runs sr where
// sr.name=LEFT(centroid,LOCATE('_',centroid) - 1));
// update clusters set member_run_id=(SELECT sr.id from sra_runs sr where
// sr.name=LEFT(`member`,LOCATE('_',`member`) - 1));

public class RNARefineryCDHitTable {

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
  private static final String QUOTE = "'";
  private static final String TAB = "\t";

  private static Logger consoleLogger;
  private static Level LOGGING_LEVEL = Level.INFO;

  static {
    consoleLogger = (Logger) LoggerFactory.getLogger(RNARefineryCDHitTable.class);
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

    bw = new BufferedWriter(new FileWriter(new File(DATA_PREFIX + "clusters_table.txt")));

    for (int i = 0; i < recs.size(); i++) {

      CDHITRec cluster = recs.get(i);

      StringBuffer sb = new StringBuffer();

      if (cluster.getEntries().isEmpty()) {
        // cluster is a singleton
        sb.append(cluster.getCentroid().getName()).append(TAB);
        sb.append("").append(TAB);
        sb.append("");
        sb.append(CR);
      }
      else {
        // cluster has members
        for (CDHITClusterEntry e : cluster.getEntries()) {
          sb.append(cluster.getCentroid().getName()).append(TAB);
          sb.append(e.getName()).append(TAB);
          sb.append(e.getSimilarity());
          sb.append(CR);
        }
      }

      bw.write(sb.toString());
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
