package net.seninp.cdhit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import org.junit.Test;

public class TestCDHITClusterEntry {

  private static final String entryString = "   3 2477nt, >ERR348576_CL353Contig1_1... at +/92.53%  ";

  @Test
  public void testParser() {
    CDHITClusterEntry ce = new CDHITClusterEntry(entryString);

    assertEquals("Testing entry length", ce.getLength(), 2477);

    assertTrue("Testing entry name", ce.getName().equals("ERR348576_CL353Contig1_1"));

    assertTrue("Testing strand", ce.getStrand().equals("+"));

    assertEquals("Testing similarity", ce.getSimilarity(), 92.53, 0.001);

  }

}
