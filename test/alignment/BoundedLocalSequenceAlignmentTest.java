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
        BoundedLocalSequenceAlignment instance = new BoundedLocalSequenceAlignment(GAP_OPEN, GAP_EXT, 1000, 5, 'N');
        StringBuilder seq1 = new StringBuilder("CTTTATGAAGAAAAAAGTT");
        StringBuilder seq2 = new StringBuilder("AAAAACTTTATGAAGAAAAAAGTTAAAAA");
        instance.align(seq1, seq2);
        int max[] = instance.get_max_coordinates();
        m = seq1.length();
        n = seq2.length();
        System.out.println("Sequence 1 of length: " + m + "\n" + seq1 + "\n");
        System.out.println("Sequence 2 of length: " + n + "\n" + seq2 + "\n");
        
        
        System.out.print("   ");
        for (j = 1; j <= 11; j++) 
            System.out.print(String.format("%3d", j ));
        System.out.println();
      //  System.out.print("   ");
//        for (j = 1; j <= 11; j++) 
  //          System.out.print(String.format("%3c", seq2.charAt(j-1) ));
    //    System.out.println();
        for (i = 1; i <= m; i++) {
            System.out.print(i + " " + seq1.charAt(i-1));
            for (j = 1; j <= 11; j++) {
                System.out.print(String.format("%3c",instance.get_direction(i,j)));
            }
            System.out.println();
        }
        System.out.println("Coordinates = "+ max[0] + " " + max[1]);
        System.out.println("Score = "+ instance.get_similarity_score());
        System.out.println(instance.get_alignment());
        instance.get_cigar();
        System.out.println(instance.get_cigar());
        System.out.println();
        
        System.out.print("   ");
        for (j = 1; j <= 11; j++) 
            System.out.print(String.format("%3d", j ));
        System.out.println();
    //    System.out.print("   ");
     //   for (j = 1; j <= n; j++) 
       //     System.out.print(String.format("%3c", seq2.charAt(j-1) ));
        //System.out.println();
        for (i = 1; i <= m; i++) {
            System.out.print(i + " " + seq1.charAt(i-1));
            for (j = 1; j <= 11; j++) {
                System.out.print(String.format("%3d", instance.get_matrix(i,j)));
            }
            System.out.println();
        }
        System.out.println("Coordinates = "+ max[0] + " " + max[1]);
        System.out.println("Score = "+ instance.get_similarity_score());
        System.out.println(instance.get_alignment());
        instance.get_cigar();
        System.out.println(instance.get_cigar());
        System.out.println();
        
        System.out.print("   ");
        for (j = 1; j <= 1; j++) 
            System.out.print(String.format("%3d", j ));
        System.out.println();
        //System.out.print("   ");
        //for (j = 1; j <= n; j++) 
        //    System.out.print(String.format("%3c", seq2.charAt(j-1) ));
        //System.out.println();
        for (i = 1; i <= m; i++) {
            System.out.print(i + " " + seq1.charAt(i-1));
            for (j = 1; j <= 11; j++) {
                System.out.print(String.format("%3d", instance.get_up(i,j)));
            }
            System.out.println();
        }        
        System.out.println("Coordinates = "+ max[0] + " " + max[1]);
        System.out.println("Score = "+ instance.get_similarity_score());
        System.out.println(instance.get_alignment());
        instance.get_cigar();
        System.out.println(instance.get_cigar());
        System.out.println();
        
        System.out.print("   ");
        for (j = 1; j <= 11; j++) 
            System.out.print(String.format("%3d", j ));
        System.out.println();
        //System.out.print("   ");
        //for (j = 1; j <= n; j++) 
        //    System.out.print(String.format("%3c", seq2.charAt(j-1) ));
        //System.out.println();
        for (i = 1; i <= m; i++) {
            System.out.print(i + " " + seq1.charAt(i-1));
            for (j = 1; j <= 11; j++) {
                System.out.print(String.format("%3d", instance.get_left(i,j)));
            }
            System.out.println();
        } 
        System.out.println("Coordinates = "+ max[0] + " " + max[1]);
        System.out.println("Score = "+ instance.get_similarity_score());
        System.out.println(instance.get_alignment());
        instance.get_cigar();
        System.out.println(instance.get_cigar());
    }
}
