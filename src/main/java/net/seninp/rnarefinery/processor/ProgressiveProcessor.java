package net.seninp.rnarefinery.processor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;
import java.util.Map.Entry;
import java.util.Set;
import net.seninp.rnarefinery.cdhit.CDHITClusterProcessor;
import net.seninp.rnarefinery.cdhit.CDHITRec;

public class ProgressiveProcessor {

  private static final DecimalFormatSymbols otherSymbols = new DecimalFormatSymbols(Locale.US);
  private static DecimalFormat df = new DecimalFormat("0.00", otherSymbols);

  private static final String prefix = "/media/Stock/RNARefinery/logs1/";

  private static final String[] INPUT_FASTA_NAMES = { "step1_t2.fa", "step2_t2.fa", "step3_t2.fa",
      "step4_t2.fa", "step5_t2.fa", "step6_t2.fa", "step7_t2.fa", "step8_t2.fa", "step9_t2.fa",
      "step10_t2.fa", "step11_t2.fa", "step12_t2.fa" };

  private static final String[] INPUT_CLUSTER_NAMES = { "step1_t2.cdhit_out.clstr",
      "step2_t2.cdhit_out.clstr", "step3_t2.cdhit_out.clstr", "step4_t2.cdhit_out.clstr",
      "step5_t2.cdhit_out.clstr", "step6_t2.cdhit_out.clstr", "step7_t2.cdhit_out.clstr",
      "step8_t2.cdhit_out.clstr", "step9_t2.cdhit_out.clstr", "step10_t2.cdhit_out.clstr",
      "step11_t2.cdhit_out.clstr", "step12_t2.cdhit_out.clstr" };

  private static final String[] SEQ_MASKS = { "ERR348560", "ERR348561", "ERR348564", "ERR348566",
      "ERR348568", "ERR348569", "ERR348571", "ERR348573", "ERR348576", "ERR348577", "ERR348579",
      "ERR348580" };

  private static final String CR = "\n";
  private static final String COMMA = ", ";
  private static final Object SPACE = " ";

  /**
   * Main runnable.
   * 
   * @param args
   * @throws Exception
   */
  public static void main(String[] args) throws Exception {

    // the data structure that keeps all tidy
    ArrayList<HashMap<String, CDHITRec>> files = new ArrayList<HashMap<String, CDHITRec>>();

    //
    // reading in all the input
    //
    {
      ArrayList<Integer> nClust = new ArrayList<Integer>();
      ArrayList<Integer> nSinglets = new ArrayList<Integer>();

      // this set takes care about non-singletons progressively
      //
      HashSet<String> notSingletons = new HashSet<String>();

      // reading the data in
      for (int k = 0; k < INPUT_CLUSTER_NAMES.length; k++) {

        String clusterFname = INPUT_CLUSTER_NAMES[k];
        System.out.println("Processing " + clusterFname);

        HashMap<String, CDHITRec> clusters = CDHITClusterProcessor.readClusters(prefix
            + clusterFname);

        nClust.add(clusters.size());
        nSinglets.add(countSingletons(clusters, notSingletons));
        files.add(clusters);
      }

      System.out.println("\n\nTop-level stats: \n");
      System.out.println("total_clusters <- "
          + CDHITClusterProcessor.arrayToRSequence(nClust.toArray(new Integer[nClust.size()])));
      System.out
          .println("total_singlets <- "
              + CDHITClusterProcessor.arrayToRSequence(nSinglets.toArray(new Integer[nSinglets
                  .size()])));

    }

    // look onto clusters centroid lengths evolution
    //
    {
      ArrayList<Double> cMeanLength = new ArrayList<Double>();
      ArrayList<Integer> cMinLength = new ArrayList<Integer>();
      ArrayList<Integer> cMaxLength = new ArrayList<Integer>();
      ArrayList<Double> sMeanLength = new ArrayList<Double>();
      ArrayList<Integer> sMinLength = new ArrayList<Integer>();
      ArrayList<Integer> sMaxLength = new ArrayList<Integer>();
      HashSet<String> notSingletons = new HashSet<String>();
      for (int i = 0; i < files.size(); i++) {
        HashMap<String, CDHITRec> entry = files.get(i);
        HashMap<String, CDHITRec> entryClusters = extractOnlyClusters(entry, notSingletons);
        saveCentroidsLength(SEQ_MASKS[i], "clusters_", entryClusters);
        HashMap<String, CDHITRec> entrySinglets = extractOnlySingletons(entry, notSingletons);
        saveCentroidsLength(SEQ_MASKS[i], "singlets_", entrySinglets);
        cMeanLength.add(CDHITClusterProcessor.meanCentroidLength(entryClusters));
        cMinLength.add(CDHITClusterProcessor.minCentroidLength(entryClusters));
        cMaxLength.add(CDHITClusterProcessor.maxCentroidLength(entryClusters));
        sMeanLength.add(CDHITClusterProcessor.meanCentroidLength(entrySinglets));
        sMinLength.add(CDHITClusterProcessor.minCentroidLength(entrySinglets));
        sMaxLength.add(CDHITClusterProcessor.maxCentroidLength(entrySinglets));
      }
      System.out.println("\n\nCentroids stats: \n");
      System.out
          .println("singletons_min_length <- "
              + CDHITClusterProcessor.arrayToRSequence(sMinLength.toArray(new Integer[sMinLength
                  .size()])));
      System.out.println("singletons_mlength <- "
          + CDHITClusterProcessor.arrayToRSequence(sMeanLength.toArray(new Double[sMeanLength
              .size()])));
      System.out
          .println("singletons_max_length <- "
              + CDHITClusterProcessor.arrayToRSequence(sMaxLength.toArray(new Integer[sMaxLength
                  .size()])));
      System.out
          .println("clusters_min_length <- "
              + CDHITClusterProcessor.arrayToRSequence(cMinLength.toArray(new Integer[cMinLength
                  .size()])));
      System.out.println("clusters_mlength <- "
          + CDHITClusterProcessor.arrayToRSequence(cMeanLength.toArray(new Double[cMeanLength
              .size()])));
      System.out
          .println("clusters_max_length <- "
              + CDHITClusterProcessor.arrayToRSequence(cMaxLength.toArray(new Integer[cMaxLength
                  .size()])));
    }

    // all the masks printed
    System.out.println("\n\nOverall counts by step: \n");
    System.out.println(Arrays.toString(SEQ_MASKS));

    // counting things by masks
    for (int i = 0; i < files.size(); i++) {
      HashMap<String, CDHITRec> entry = files.get(i);
      StringBuffer sb = new StringBuffer();
      HashMap<String, Integer> counts = countAllByMask(entry);
      for (String mask : SEQ_MASKS) {
        Integer count = counts.get(mask);
        if (null == count) {
          sb.append("0, ");
        }
        else {
          sb.append(count).append(COMMA);
        }
      }
      System.out.println(sb.toString());
    }

    // look onto singletons evolution
    //
    System.out.println("\n\nSingletons evolution: \n");
    {
      ArrayList<HashMap<String, CDHITRec>> singletons = new ArrayList<HashMap<String, CDHITRec>>();
      HashSet<String> notSingletons = new HashSet<String>();
      // we keep the following data structure to mark non-singleton sequences
      for (int i = 0; i < files.size(); i++) {
        HashMap<String, CDHITRec> entry = files.get(i);
        HashMap<String, CDHITRec> entrySingletons = extractOnlySingletons(entry, notSingletons);
        singletons.add(entrySingletons);
      }
      //
      //
      System.out.println(Arrays.toString(SEQ_MASKS));
      for (int i = 0; i < singletons.size(); i++) {
        String mask = SEQ_MASKS[i];
        HashMap<String, CDHITRec> startingSet = singletonsByMask(singletons.get(i), mask);
        for (int j = 0; j < i; j++) {
          System.out.print("0, ");
        }
        System.out.print(startingSet.size() + ", ");
        for (int k = i + 1; k < singletons.size(); k++) {
          HashMap<String, CDHITRec> cSet = singletons.get(k);
          // intersect and get result
          System.out.print(intersectAndCount(startingSet.keySet(), cSet.keySet()) + ",");
        }
        System.out.println();
      }
    }

    // look onto clusters evolution
    //
    System.out.println("\n\nClusters evolution: \n");
    {
      HashSet<String> notSingletons = new HashSet<String>();
      ArrayList<HashMap<String, CDHITRec>> clusters = new ArrayList<HashMap<String, CDHITRec>>();
      for (int i = 0; i < files.size(); i++) {
        HashMap<String, CDHITRec> entry = files.get(i);
        HashMap<String, CDHITRec> entrySingletons = extractOnlyClusters(entry, notSingletons);
        clusters.add(entrySingletons);
      }
      //
      //
      System.out.println(Arrays.toString(SEQ_MASKS));
      for (int i = 0; i < clusters.size(); i++) {
        String mask = SEQ_MASKS[i];
        HashMap<String, CDHITRec> startingSet = clustersByMask(clusters.get(i), mask);
        for (int j = 0; j < i; j++) {
          System.out.print("0, ");
        }
        System.out.print(startingSet.size() + ", ");
        for (int k = i + 1; k < clusters.size(); k++) {
          HashMap<String, CDHITRec> cSet = clusters.get(k);
          // intersect and get result
          System.out.print(intersectAndCount(startingSet.keySet(), cSet.keySet()) + ",");
        }
        System.out.println();
      }
    }

  }

  private static void saveCentroidsLength(String mask, String prefix,
      HashMap<String, CDHITRec> entry) throws IOException {
    BufferedWriter bw = new BufferedWriter(new FileWriter(new File(prefix + mask + ".txt")));
    for (Entry<String, CDHITRec> e : entry.entrySet()) {
      bw.write(e.getValue().getCentroid().getLength() + "\n");
    }
    bw.close();
  }

  /**
   * Get only clusters.
   * 
   * @param entry
   * @param notSingletons
   * @return
   */
  private static HashMap<String, CDHITRec> extractOnlyClusters(HashMap<String, CDHITRec> entry,
      HashSet<String> notSingletons) {
    HashMap<String, CDHITRec> res = new HashMap<String, CDHITRec>();
    for (Entry<String, CDHITRec> e : entry.entrySet()) {
      if (1 != e.getValue().getSize()) {
        res.put(e.getKey(), e.getValue());
        notSingletons.add(e.getKey());
      }
      else if (notSingletons.contains(e.getKey())) {
        res.put(e.getKey(), e.getValue());
      }
    }
    return res;
  }

  /**
   * Extracts only singletons out of the CDHit collection.
   * 
   * @param entry The collection.
   * @param notSingletons
   * @return only singletons sub-collection.
   */
  private static HashMap<String, CDHITRec> extractOnlySingletons(HashMap<String, CDHITRec> entry,
      HashSet<String> notSingletons) {
    HashMap<String, CDHITRec> res = new HashMap<String, CDHITRec>();
    for (Entry<String, CDHITRec> e : entry.entrySet()) {
      String cName = e.getKey();
      if ((1 == e.getValue().getSize()) && (!notSingletons.contains(cName))) {
        res.put(e.getKey(), e.getValue());
      }
      else {
        notSingletons.add(e.getKey());
      }
    }
    return res;
  }

  /**
   * Get only clusters that match the mask.
   * 
   * @param entry
   * @param mask
   * @return
   */
  private static HashMap<String, CDHITRec> clustersByMask(HashMap<String, CDHITRec> entry,
      String mask) {
    HashMap<String, CDHITRec> res = new HashMap<String, CDHITRec>();
    for (Entry<String, CDHITRec> e : entry.entrySet()) {
      if (mask.equalsIgnoreCase(e.getKey().substring(0, 9))) {
        res.put(e.getKey(), e.getValue());
      }
    }
    return res;
  }

  /**
   * Get only singletons that match the mask.
   * 
   * @param entry
   * @param mask
   * @return
   */
  private static HashMap<String, CDHITRec> singletonsByMask(HashMap<String, CDHITRec> entry,
      String mask) {
    HashMap<String, CDHITRec> res = new HashMap<String, CDHITRec>();
    for (Entry<String, CDHITRec> e : entry.entrySet()) {
      if (mask.equalsIgnoreCase(e.getKey().substring(0, 9))) {
        res.put(e.getKey(), e.getValue());
      }
    }
    return res;
  }

  /**
   * Counts singletons in a collection.
   * 
   * @param clusters the clusters collection.
   * @param notSingletons
   * @return counts.
   */
  private static int countSingletons(HashMap<String, CDHITRec> clusters,
      HashSet<String> notSingletons) {
    int res = 0;
    for (Entry<String, CDHITRec> e : clusters.entrySet()) {
      String cName = e.getKey();
      if ((1 == e.getValue().getSize()) && (!notSingletons.contains(cName))) {
        res++;
      }
      else {
        notSingletons.add(cName);
      }
    }
    return res;
  }

  /**
   * Counts centroids by the dataset mask.
   * 
   * @param the clusters map.
   * @return counts.
   */
  private static HashMap<String, Integer> countAllByMask(HashMap<String, CDHITRec> entry) {
    HashMap<String, Integer> res = new HashMap<String, Integer>();
    for (Entry<String, CDHITRec> e : entry.entrySet()) {
      String centroidMask = e.getKey().substring(0, 9);
      if (!(res.containsKey(centroidMask))) {
        res.put(centroidMask, 1);
      }
      else {
        res.put(centroidMask, res.get(centroidMask) + 1);
      }
    }
    return res;
  }

  /**
   * Intersects two sets.
   * 
   * @param keySet
   * @param keySet2
   * @return
   */
  private static int intersectAndCount(Set<String> keySet, Set<String> keySet2) {
    Set<String> intersection = new HashSet<String>(keySet); // use the copy constructor
    intersection.retainAll(new HashSet<String>(keySet2));
    return intersection.size();
  }

}
