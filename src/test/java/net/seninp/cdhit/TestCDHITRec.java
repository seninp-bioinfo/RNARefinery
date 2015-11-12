package net.seninp.cdhit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import org.junit.Test;

public class TestCDHITRec {

  private static final String buffer = ">Cluster 0\n"
      + "0 22838nt, >ERR348569_k25_Locus_171_Trans... *\n"
      + "1 967nt, >ERR348576_CL11246Contig1_1... at -/100.00%\n"
      + "2 737nt, >ERR348576_CL14788Contig1_1... at +/99.86%\n"
      + "3 2477nt, >ERR348576_CL353Contig1_1... at +/92.53%\n"
      + "4 1734nt, >ERR348576_k25_Locus_32677_Tra... at -/100.00%\n"
      + "5 346nt, >ERR348576_k25_Locus_48053_Tra... at +/100.00%\n"
      + "6 1437nt, >ERR348576_k43_Locus_51503_Tra... at -/100.00%\n";

  @Test
  public void testParser() {
    CDHITRec rec = CDHITRec.fromString(buffer);
    assertEquals("testing id", 0, rec.getId());

    CDHITClusterEntry centroid = rec.getCentroid();
    assertEquals("testing centroid", 22838, centroid.getLength());
    assertEquals("testing centroid", 100.0, centroid.getSimilarity(), 0.001);
    assertTrue("testing centroid", "+".equals(centroid.getStrand()));
    assertTrue("testing centroid", "ERR348569_k25_Locus_171_Trans".equals(centroid.getName()));
  }

}
