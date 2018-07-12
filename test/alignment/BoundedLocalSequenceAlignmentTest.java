/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package alignment;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static pantools.Pantools.GAP_EXT;
import static pantools.Pantools.GAP_OPEN;

/**
 *
 * @author sheik005
 */
public class BoundedLocalSequenceAlignmentTest {
    
    public BoundedLocalSequenceAlignmentTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }

    /**
     * Test of initialize_NUCC_matrix method, of class LocalSequenceAlignment.
     */
    @Test
    public void test() {
        int i, j, m, n;
        BoundedLocalSequenceAlignment instance = new BoundedLocalSequenceAlignment(GAP_OPEN, GAP_EXT, 1000, 3, true, 'N');
        instance.initialize_bound(3, 1000);
        StringBuilder seq1, seq2;
        seq2 = new StringBuilder("GGGAAAAACTTTATGAAGAAAAAAGTTAAAAAGGG");
        seq1 = new StringBuilder(   "AAAAACTTTATGAAGAAAAAAGTTAAAAA");
        instance.align(seq1, seq2);
        seq2 = new StringBuilder("GGGGGGAAAAACTTTATGAAGAAAAAAGTTAAAAA");
        seq1 = new StringBuilder(   "AAAAACTTTATGAAGAAAAAAGTTAAAAA");
        instance.align(seq1, seq2);
        seq2 = new StringBuilder("GGGGAAAAACTTTATGAAGAAAAAAGTTAAAAAGG");
        seq1 = new StringBuilder(   "AAAAACTTTATGAAGAAAAAAGTTAAAAA");
        instance.align(seq1, seq2);
        seq2 = new StringBuilder("GGGGAAAAACTTTATGAAGAAAGTTAAAAAGGGG");
        seq1 = new StringBuilder(   "AAAAACTTTATGAAGAAAAAAGTTAAAA");
        instance.align(seq1, seq2);
        seq2 = new StringBuilder("GGGAAAAACTTTATGAAGAAAAAAGTTAAAAA");
        seq1 = new StringBuilder(   "AAAAACTTTATGAAGAAAGTTAAAAA");
        instance.align(seq1, seq2);
        seq2 = new StringBuilder("GGGAAAAACTTTATGAAGAAAAAAGTTAAAAAGGGGG");
        seq1 = new StringBuilder(   "GGAAAAACTTTATGAAGAAAAAAGTTAAAAA");
        instance.align(seq1, seq2);
        seq2 = new StringBuilder("GGGAAAAACTTTATGAAGAAAAAAGTTAAAAAGGGGG");
        seq1 = new StringBuilder(   "AAAAACTTTATGAAGAAAAAAGTTAAAAAGG");
        instance.align(seq1, seq2);    }
}
