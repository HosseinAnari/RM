/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pangenome;

import alignment.BoundedLocalSequenceAlignment;
import alignment.LocalSequenceAlignment;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import sequence.SequenceDatabase;
import sequence.SequenceScanner;
import index.IndexPointer;
import index.IndexDatabase;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Relationship;
import org.neo4j.graphdb.RelationshipType;
import org.neo4j.graphdb.ResourceIterator;
import org.neo4j.graphdb.Transaction;
import org.neo4j.graphdb.factory.GraphDatabaseFactory;
import org.neo4j.graphdb.Direction;
import static org.neo4j.graphdb.factory.GraphDatabaseSettings.keep_logical_logs;
import org.neo4j.io.fs.FileUtils;
import static pantools.Pantools.GENOME_DATABASE_PATH;
import static pantools.Pantools.GRAPH_DATABASE_PATH;
import static pantools.Pantools.INDEX_DATABASE_PATH;
import pantools.Pantools.RelTypes;
import static pantools.Pantools.PATH_TO_THE_PANGENOME_DATABASE;
import static pantools.Pantools.PATH_TO_THE_GENOMES_FILE;
import static pantools.Pantools.PATH_TO_THE_REGIONS_FILE;
import static pantools.Pantools.PATH_TO_THE_GENOME_NUMBERS_FILE;
import static pantools.Pantools.RAW_ABUNDANCE_FILE;
import static pantools.Pantools.num_bases;
import static pantools.Pantools.num_degenerates;
import static pantools.Pantools.num_edges;
import static pantools.Pantools.num_nodes;
import static pantools.Pantools.phaseTime;
import static pantools.Pantools.startTime;
import static pantools.Pantools.MAX_TRANSACTION_SIZE;
import static pantools.Pantools.ANCHORS_DISTANCE;
import static pantools.Pantools.complement;
import static pantools.Pantools.write_fasta;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Date;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.PriorityQueue;
import java.util.Queue;
import java.util.Random;
import java.util.Scanner;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import org.neo4j.graphdb.NotFoundException;
import static pantools.Pantools.DEBUG;
import static pantools.Pantools.GAP_EXT;
import static pantools.Pantools.GAP_OPEN;
import static pantools.Pantools.K_SIZE;
import static pantools.Pantools.PATH_TO_THE_FIRST_SRA;
import static pantools.Pantools.PATH_TO_THE_SECOND_SRA;
import static pantools.Pantools.SHOW_KMERS;
import static pantools.Pantools.CLIPPING_STRINGENCY;
import static pantools.Pantools.BAMFORMAT;
import static pantools.Pantools.OUTPUT_PATH;
import static pantools.Pantools.THREADS;
import static pantools.Pantools.degenerate_label;
import static pantools.Pantools.genome_label;
import static pantools.Pantools.low_complexity_label;
import static pantools.Pantools.nucleotide_label;
import static pantools.Pantools.pangenome_label;
import static pantools.Pantools.sequence_label;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.logging.Level;
import java.util.logging.Logger;
import static pantools.Pantools.ALIGNMENT_BOUND;
import static pantools.Pantools.ALIGNMENT_MODE;
import static pantools.Pantools.NUM_KMER_SAMPLES;
import static pantools.Pantools.MAX_ALIGNMENT_LENGTH;
import static pantools.Pantools.MAX_FRAGMENT_LENGTH;
import static pantools.Pantools.MIN_HIT_LENGTH;
import static pantools.Pantools.MIN_IDENTITY;
import static pantools.Pantools.MAX_NUM_LOCATIONS;
import static pantools.Pantools.genomeDb;
import static pantools.Pantools.genome_scanner;
import static pantools.Pantools.graphDb;
import static pantools.Pantools.indexDb;

/**
 * Implements all the functionalities related to the sequence layer of the pangenome
 * 
 * @author Siavash Sheikhizadeh, Bioinformatics chairgroup, Wageningen
 * University, Netherlands
 */
public class GenomeLayer {

    private Node curr_node;
    private byte curr_side;        
    private boolean finish;
    private AtomicInteger[] num_shared_mapping;
    private AtomicInteger[] num_unique_mapping;
    private AtomicInteger[] num_unmapped;
    private AtomicInteger number_of_alignments;
    private AtomicInteger number_of_hits;
    private AtomicInteger read_number;
    private SequenceScanner genome_scanner;    
    
    public class read{
        StringBuilder name;
        StringBuilder forward_seq;
        StringBuilder reverse_seq;
        //StringBuilder quality;
        
        public read(){
            name = new StringBuilder();
            forward_seq = new StringBuilder();
            reverse_seq = new StringBuilder();
            //quality = new StringBuilder();
        }
        
        public void clear(){
            name.setLength(0);
            forward_seq.setLength(0);
            reverse_seq.setLength(0);
            //quality.setLength(0);
        }
        
        public int length(){
            return forward_seq.length();
        }
    }
    
    public class single_hit{
        public int genome;
        public int sequence;
        public double identity;
        public int start;
        public int offset;
        public int length;
        public int deletions;
        public boolean forward;
        public String cigar;
        public String reference;
        public single_hit(int gn, int sq, double idn, int st, int off, int len, int del, boolean fw, String cg, String r){
            genome = gn;
            sequence = sq;
            identity = idn;
            start = st;
            offset = off;
            length = len;
            deletions = del;
            forward = fw;
            cigar = cg;
            reference = r;
        }
        
        public single_hit(single_hit h){
            genome = h.genome;
            sequence = h.sequence;
            identity = h.identity;
            start = h.start;
            offset = h.offset;
            length = h.length;
            deletions = h.deletions;
            forward = h.forward;
            cigar = h.cigar;
            reference = h.reference;
        }
        

        @Override
        public String toString(){
            return "(genome:" + genome + 
                    ",sequence:" + sequence + 
                    ",identity:" + identity + 
                    ",start:" + start + 
                    ",offset:" + offset + 
                    ",length:" + length + 
                    ",deletions:" + deletions + 
                    ",forward:" + forward + 
                    ",reference:" + reference + 
                    ",cigar:" + cigar +")";
        }
    }

    public class paired_hit{
        public int fragment_length;
        
        public int genome1;
        public int sequence1;
        public double identity1;
        public int start1;
        public int offset1;
        public int length1;
        public int deletions1;
        public boolean forward1;
        public String cigar1;
        public String reference1;

        public int genome2;
        public int sequence2;
        public double identity2;
        public int start2;
        public int offset2;
        public int length2;
        public int deletions2;
        public boolean forward2;
        public String cigar2;
        public String reference2;
        
        public paired_hit(int flen, int gn1, int sq1, double idn1, int st1, int off1, int len1, int del1, boolean fw1, String cg1, String r1,
                          int gn2, int sq2, double idn2, int st2, int off2, int len2, int del2, boolean fw2, String cg2, String r2){
            fragment_length = flen;
            genome1 = gn1;
            sequence1 = sq1;
            identity1 = idn1;
            start1 = st1;
            offset1 = off1;
            length1 = len1;
            deletions1 = del1;
            forward1 = fw1;
            cigar1 = cg1;
            reference1 = r1;

            genome2 = gn2;
            sequence2 = sq2;
            identity2 = idn2;
            start2 = st2;
            offset2 = off2;
            length2 = len2;
            deletions2 = del2;
            forward2 = fw2;
            cigar2 = cg2;
            reference2 = r2;
        }

        public paired_hit(paired_hit h){
            fragment_length = h.fragment_length;
            genome1 = h.genome1;
            sequence1 = h.sequence1;
            identity1 = h.identity1;
            start1 = h.start1;
            offset1 = h.offset1;
            length1 = h.length1;
            deletions1 = h.deletions1;
            forward1 = h.forward1;
            cigar1 = h.cigar1;
            reference1 = h.reference1;

            genome2 = h.genome2;
            sequence2 = h.sequence2;
            identity2 = h.identity2;
            start2 = h.start2;
            offset2 = h.offset2;
            length2 = h.length2;
            deletions2 = h.deletions2;
            forward2 = h.forward2;
            cigar2 = h.cigar2;
            reference2 = h.reference2;
        }
        

        public double get_identity(){
            return identity1 + identity2;
        }

        public int get_min_start(){
            return Math.min(start1, start2);
        }

        public int get_max_start(){
            return Math.max(start1, start2);
        }

        @Override
        public String toString(){
            return "(genome1:" + genome1 + 
                    ",sequence1:" + sequence1 + 
                    ",identity1:" + identity1 + 
                    ",start1:" + start1 + 
                    ",offset1:" + offset1 + 
                    ",length1:" + length1 + 
                    ",deletions1:" + deletions1 + 
                    ",forward1:" + forward1 + 
                    ",reference1:" + reference1 + 
                    ",cigar1:" + cigar1 +")" +
                    "\n" +
                    "(genome2:" + genome2 + 
                    ",sequence2:" + sequence2 + 
                    ",identity2:" + identity2 + 
                    ",start2:" + start2 + 
                    ",offset2:" + offset2 + 
                    ",length2:" + length2 + 
                    ",deletions2:" + deletions2 + 
                    ",forward2:" + forward2 + 
                    ",reference2:" + reference2 + 
                    ",cigar2:" + cigar2 +")";
        }
    }

    /**
     * Implements a comparator for integer arrays of size two
     */
    public static class single_hitComparator implements Comparator<single_hit> {
        @Override
        public int compare(single_hit x, single_hit y) {
            if (y.identity < x.identity) 
                return -1;
            else if (y.identity > x.identity) 
                return 1;
            else if (y.cigar != null && x.cigar != null){
                 if (y.cigar.length() > x.cigar.length()) 
                    return -1;
                else if (y.cigar.length() < x.cigar.length()) 
                    return 1;
            } 
            if (y.genome > x.genome) 
                return -1;
            else if (y.genome < x.genome) 
                return 1;
            else if (y.sequence > x.sequence) 
                return -1;
            else if (y.sequence < x.sequence) 
                return 1;
            else if (y.start < x.start) 
                return -1;
            else if (y.start > x.start) 
                return 1;
            else
                return 0;
        }
    }      

    public static class paired_hitComparator implements Comparator<paired_hit> {
        @Override
        public int compare(paired_hit x, paired_hit y) {
            if (y.get_identity() < x.get_identity()) 
                return -1;
            else if (y.get_identity() > x.get_identity()) 
                return 1;
            else if (y.fragment_length > x.fragment_length) 
                return -1;
            else if (y.fragment_length < x.fragment_length) 
                return 1;
            else if (y.cigar1 != null && x.cigar1 != null && y.cigar2 != null && x.cigar2 != null){
                if (y.cigar1.length() + y.cigar2.length() > x.cigar1.length() + x.cigar2.length()) 
                    return -1;
                else if (y.cigar1.length() + y.cigar2.length() < x.cigar1.length() + x.cigar2.length()) 
                    return 1;
            } 
            if (y.genome1 > x.genome1) 
                return -1;
            else if (y.genome1 < x.genome1) 
                return 1;
            else if (y.sequence1 > x.sequence1) 
                return -1;
            else if (y.sequence1 < x.sequence1) 
                return 1;
            else if (y.start1 < x.start1) 
                return -1;
            else if (y.start1 > x.start1) 
                return 1;
            else
                return 0;
        }
    }      

    
    public static class CountComparator implements Comparator<int[]> {
        @Override
        public int compare(int[] x, int[] y) {
            if (y[1] < x[1]) 
                return -1;
            else if (y[1] > x[1]) 
                return 1;
            else
                return 0;
        }
    }    

    public static class IntComparator implements Comparator<Integer> {

        @Override
        public int compare(Integer v1, Integer v2) {
            return v1 < v2 ? -1 : v1 > v2 ? +1 : 0;
        }
    }

    public class Map implements Runnable {
        SequenceDatabase genome_database;
        IndexDatabase index_database;
        int thread_id;
        int K;
        read[] reads;
        SequenceScanner[] reads_scanner;    
        ArrayList<Integer>[][][] locations;
        IntComparator intcomp = new IntComparator();
        IndexPointer anchor_pointer;
        StringBuilder reference;
        LinkedList<Integer>[][] sequences;
        BoundedLocalSequenceAlignment bounded_aligner;
        LocalSequenceAlignment aligner;
        PriorityQueue<single_hit>[][] hits;
        LinkedList<single_hit>[][] alignments;
        PriorityQueue<single_hit> single_hits;
        Queue<single_hit> single_hits_2;
        PriorityQueue<paired_hit> paired_hits;
        Queue<paired_hit> paired_hits_2;
        int anchor_position;
        boolean paired_end;
        CountComparator cnt_comp;
        ArrayList<int[]> hit_counts;
        int num_neighbors = 0;
        int num_kmers = 0;
        int num_found = 0;
        int[] address = new int[3];
        int num_segments;
        LinkedList<Integer> genome_numbers;
        SAMFileWriter[] sam_writers;
        SAMFileHeader[] sam_headers;
        SequenceDatabase[] sequencingDb;
        long[] genome_sizes;
        int[] shared;
        int[] unique;
        int[] unmapped;
        int total_unique = 1;
        int num_genomes;
        int[] num_sequences;
        long[][] sequence_length;
        String[][] sequence_titles;
        Random rand;
        double[] raw_abundance;
        
        public Map(int id, LinkedList gn, boolean paired, SAMFileWriter[] sams, SAMFileHeader[] headers, SequenceDatabase[] sd) {
            int i, j, genome, abun;
            index_database = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH, "sorted");
            genome_database = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH);
            sequencingDb = sd;
            K = index_database.get_K();
            num_genomes = genome_database.num_genomes;
            num_sequences = new int[num_genomes + 1];
            sequence_length = new long[num_genomes + 1][];
            sequence_titles = new String[num_genomes + 1][];
            genome_scanner = new SequenceScanner(genome_database, 1, 1, K, index_database.get_pre_len());
            sam_writers = sams;
            sam_headers = headers;
            thread_id = id;
            genome_numbers = gn;
            paired_end = paired; 
            hits = new PriorityQueue[2][];
            alignments = new LinkedList[2][];
            reads_scanner = new SequenceScanner[2];
            sequences = new LinkedList[2][];
            locations = new ArrayList[2][][];
            reads = new read[2];
            
            hits[0] = new PriorityQueue[num_genomes + 1];
            alignments[0] = new LinkedList[num_genomes + 1];
            single_hits = new PriorityQueue(new single_hitComparator());
            single_hits_2 = new LinkedList();
            paired_hits = new PriorityQueue(new paired_hitComparator());
            paired_hits_2 = new LinkedList();
            reads_scanner[0] = new SequenceScanner(sequencingDb[0], 1, 1, K, index_database.get_pre_len());
            sequences[0] = new LinkedList[num_genomes + 1];
            locations[0] = new ArrayList[num_genomes + 1][];
            genome_sizes = new long[num_genomes + 1];
            shared = new int[num_genomes + 1];
            unique = new int[num_genomes + 1];
            unmapped = new int[num_genomes + 1];
            raw_abundance = new double[num_genomes + 1];
            reads[0] = new read();
            for (i = 0; i <= num_genomes; ++ i){
                hits[0][i] = null;
                alignments[0][i] = null;
                locations[0][i] = null;
                sequences[0][i] = null;
                genome_sizes[i] = genome_database.genome_length[i];
                num_sequences[i] = genome_database.num_sequences[i];
                sequence_length[i] = new long[num_sequences[i] + 1];
                sequence_titles[i] = new String[num_sequences[i] + 1];
                shared[i] = 0;
            }
            for (ListIterator<Integer> genome_itr = genome_numbers.listIterator();genome_itr.hasNext();){
                genome = genome_itr.next();
                hits[0][genome] = new PriorityQueue(new single_hitComparator());
                alignments[0][genome] = new LinkedList();
                locations[0][genome] = new ArrayList[num_sequences[genome] + 1];
                for (j = 1; j <= num_sequences[genome]; ++j){
                    locations[0][genome][j] = new ArrayList();
                    sequence_length[genome][j] = genome_database.sequence_length[genome][j];
                    sequence_titles[genome][j] = genome_database.sequence_titles[genome][j];
                }
                sequences[0][genome] = new LinkedList();
            }
            if (paired_end){
                hits[1] = new PriorityQueue[num_genomes + 1];
                alignments[1] = new LinkedList[num_genomes + 1];
                reads_scanner[1] = new SequenceScanner(sequencingDb[1], 1, 1 ,K, index_database.get_pre_len());
                sequences[1] = new LinkedList[num_genomes + 1];
                locations[1] = new ArrayList[num_genomes + 1][];
                reads[1] = new read();
                for (i = 0; i <= num_genomes; ++ i){
                    hits[1][i] = null;
                    alignments[1][i] = null;
                    locations[1][i] = null;
                    sequences[1][i] = null;
                }
                for (ListIterator<Integer> genome_itr = genome_numbers.listIterator();genome_itr.hasNext();){
                    genome = genome_itr.next();
                    hits[1][genome] = new PriorityQueue(new single_hitComparator());
                    alignments[1][genome] = new LinkedList();
                    locations[1][genome] = new ArrayList[num_sequences[genome] + 1];
                    for (j = 1; j <= num_sequences[genome]; ++j)
                        locations[1][genome][j] = new ArrayList();
                    sequences[1][genome] = new LinkedList();
                }
            }
            anchor_pointer = new IndexPointer();
            reference = new StringBuilder();
            //bounded_aligner = new LocalSequenceAlignment(GAP_OPEN, GAP_EXT, MAX_ALIGNMENT_LENGTH, CLIP, 'N');
            bounded_aligner = new BoundedLocalSequenceAlignment(GAP_OPEN, GAP_EXT, MAX_ALIGNMENT_LENGTH, ALIGNMENT_BOUND, CLIPPING_STRINGENCY, 'N');
            aligner = new LocalSequenceAlignment(GAP_OPEN, GAP_EXT, MAX_ALIGNMENT_LENGTH, CLIPPING_STRINGENCY, 'N');
            cnt_comp = new CountComparator();
            hit_counts = new ArrayList();
            rand = new Random();
            if (RAW_ABUNDANCE_FILE.equals("")){
                for (i = 1; i < num_genomes; ++i)
                    raw_abundance[i] = 1.0;
            } else {
                try{
                    String line;
                    String[] fields;
                    BufferedReader in = new BufferedReader(new FileReader(RAW_ABUNDANCE_FILE));
                    while (in.ready()){
                        line = in.readLine().trim();
                        if (line.equals("") || line.startsWith("Genome"))
                        continue;
                        fields = line.split("\\s+");
                        abun = Integer.parseInt(fields[2]);
                        if (abun <= 0)
                            abun = 1;
                        raw_abundance[Integer.parseInt(fields[0])] = abun;
                    }
                    in.close();
                } catch (IOException ex){
                    
                }
            }
        }

        @Override
        public void run() {
            int read_num, trsc;
            int i;
            int num_reads;
            int counter;
            num_reads = sequencingDb[0].num_sequences[1];
            counter = thread_id * num_reads / 50 / THREADS;
            read_num = read_number.getAndIncrement();
            num_segments = paired_end ? 2 : 1;
            while (read_num <= num_reads) {
                try(Transaction tx = graphDb.beginTx()){
                    for (trsc = 1; trsc <= 10000 && read_num <= num_reads; ++trsc, ++counter, read_num = read_number.getAndIncrement()){
                        find_locations(read_num);
                        cluster_and_examine_locations();
                        report_hits();
                        if (counter % (1 + num_reads / 50) == 0)
                            //System.out.println(unique[1] + " " + unique[2]);
                            System.out.print("|");
                    }
                    tx.success();
                }
            }// for each read
            for (ListIterator<Integer> genome_itr = genome_numbers.listIterator();genome_itr.hasNext();){
                i = genome_itr.next();
                num_shared_mapping[i].getAndAdd(shared[i]);
                num_unique_mapping[i].getAndAdd(unique[i]);
                num_unmapped[i].getAndAdd(unmapped[i]);
            }
        }
        
        public void get_read(int read_num, int mate){
            int i;
            char ch;
            reads[mate].clear();
            reads_scanner[mate].set_sequence(read_num);
            reads_scanner[mate].get_complete_sequence(reads[mate].forward_seq, 1, read_num, true);
            reads_scanner[mate].get_complete_sequence(reads[mate].reverse_seq, 1, read_num, false);
            //System.out.println(reads[mate].forward_seq.toString());
            //reads_scanner[mate].get_sequence_quality(reads[mate].quality, 1, read_num);
            reads_scanner[mate].get_sequence_title(reads[mate].name, 1, read_num);
            //read_name = reads[mate].name;
            ch = reads[mate].name.charAt(reads[mate].name.length() - 1);
            if (paired_end && (ch == '1' || ch == '2'))
                reads[mate].name.setLength(reads[mate].name.length() - 2);
            if (reads[mate].length() > MAX_ALIGNMENT_LENGTH){
                System.out.println(reads[mate].length() + " exceeds the MAX_ALIGNMENT_LENGTH " + MAX_ALIGNMENT_LENGTH +
                        ". Try a larger MAX_ALIGNMENT_LENGTH with --max-length option.");
                System.exit(1);
                
            }
            if ((i = reads[mate].name.indexOf(" ")) != -1)
                reads[mate].name.setLength(i);
            if (DEBUG) 
                System.out.println("read: " + reads[mate].name);         
        }

        public void find_locations(int read_num){
            Node node = null;
            int i, mate, num_segments, step;
            long cur_index, prev_node_id;
            num_segments = paired_end ? 2 : 1;
            for (mate = 0; mate < num_segments; ++mate){
                get_read(read_num, mate);
                anchor_position = reads_scanner[mate].initialize_kmer(0);
                prev_node_id = 0;
                step = reads[mate].length() / NUM_KMER_SAMPLES;
                step = (step == 0 ? 1 : step);
                while (anchor_position != Integer.MIN_VALUE && anchor_position < reads[mate].length() - 11){
                    cur_index = reads_scanner[mate].find_curr_kmer(index_database);
                    try{
                        if (cur_index != -1l){
                            index_database.get_pointer(anchor_pointer, cur_index);
                            if (prev_node_id != anchor_pointer.node_id){ // to save a bit of time
                                node = graphDb.getNodeById(anchor_pointer.node_id);
                                prev_node_id = anchor_pointer.node_id;
                            } 
                            explore_node(node, mate);
                        }
                        for (i = 0; i < step; ++i)
                            anchor_position = reads_scanner[mate].next_kmer();
                    } catch (NotFoundException|ClassCastException ex){
                        //num_exceptions++;
                        //System.out.println(ex.getMessage());
                    } 
                }
            }
        }
        
        public void explore_node(Node node, int mate) {
            int i, loc, offset, read_len = reads[mate].length();
            long seq_len;
            char side;
            int[] location_array;
            int genome, sequence;
            boolean canonical_kmer;
            offset = anchor_pointer.offset;
            canonical_kmer = reads_scanner[mate].get_curr_kmer().get_canonical();
            int node_len;
        // for each incoming edge to the node of the anchor    
            for (Relationship r: node.getRelationships(Direction.INCOMING)){
                //num_neighbors++;
                side = r.getType().name().charAt(1);
            // for all seuences passing that node     
                for (String seq_id: r.getPropertyKeys()){
                    extract_address(address, seq_id);
                    genome = address[0];
                    if (locations[mate][genome] != null){// should map against this genome
                        sequence = address[1];
                        if (!sequences[mate][genome].contains(sequence))
                            sequences[mate][genome].add(sequence);
                    // calculate the locations based on the offsets in the node    
                        location_array = (int[])r.getProperty(seq_id);
                        seq_len = sequence_length[genome][sequence];
                        if (side == 'F'){
                            for (i = 0; i <= location_array.length - 1; i += 1){
                                if (anchor_pointer.canonical ^ canonical_kmer){
                                    loc = location_array[i] + offset - read_len + anchor_position + K;
                                    if (loc >= 0 && loc <= seq_len - read_len){
                                        locations[mate][genome][sequence].add( -(1 + loc) );
                                    }
                                    //System.out.println("F-" + loc);
                                } else {
                                    loc = location_array[i] + offset - anchor_position;
                                    if (loc >= 0 && loc <= seq_len - read_len){
                                        locations[mate][genome][sequence].add(loc);
                                    }
                                    //System.out.println("F+" + loc);
                                }
                            }
                        }else{
                            node_len = (int)node.getProperty("length");
                            for (i = 0; i <= location_array.length - 1; i += 1){
                                if (anchor_pointer.canonical ^ canonical_kmer){
                                    loc = location_array[i] + node_len - K - offset - anchor_position;
                                    if (loc >= 0 && loc <= seq_len - read_len){
                                        locations[mate][genome][sequence].add(loc);
                                    }
                                    //System.out.println("R+" + loc);
                                } else {
                                    loc = location_array[i] + node_len - offset - read_len + anchor_position;
                                    if (loc >= 0 && loc <= seq_len - read_len){
                                        locations[mate][genome][sequence].add( -(1 + loc) );
                                    }
                                    //System.out.println("R-" + loc);
                                }
                            }
                        }
                    }
                }
            }            
        }
        
        public void extract_address(int[] a, String prp){
            int i;
            char ch;
            a[0] = 0;
            a[1] = 0;
            for (i = 1; i < prp.length(); ++i)
                if ((ch = prp.charAt(i)) != 'S')
                    a[0] = a[0] * 10 + ch - 48;
                else
                    break;
            for (++i; i < prp.length(); ++i)
                a[1] = a[1] * 10 + prp.charAt(i) - 48;
        }

        public void cluster_and_examine_locations(){
            int genome, sequence, mate, first_start, prev_start;
            int ref_start, i, n, count;
            for (mate = 0; mate < num_segments; ++mate){
                for (ListIterator<Integer> genome_itr = genome_numbers.listIterator();genome_itr.hasNext();){
                    genome = genome_itr.next();
                    while (!sequences[mate][genome].isEmpty()){
                        hit_counts.clear();
                        sequence = sequences[mate][genome].remove();
                        if (locations[mate][genome][sequence].size() > 0){
                            locations[mate][genome][sequence].sort(intcomp);
                            first_start = prev_start = locations[mate][genome][sequence].get(0);
                            n = locations[mate][genome][sequence].size();
                            for (count = 0, i = 0; i < n; ++i) {
                                ref_start = locations[mate][genome][sequence].get(i);
                                //System.out.print(ref_start + " ");
                                if (ref_start - prev_start != 0){ 
                                    hit_counts.add(new int[]{first_start, count});
                                    count = 1;
                                    first_start = ref_start;
                                } else
                                    ++count;
                                prev_start = ref_start;
                            }
                            //System.out.println();
                            hit_counts.add(new int[]{first_start, count});
                            hit_counts.sort(cnt_comp);
                            n = Math.min(hit_counts.size(), MAX_NUM_LOCATIONS);//hit_counts.get(0)[1] / 2;
                            //n = hit_counts.size();
                            for (i = 0; i < n; ++i)
                                examine_location(mate, genome, sequence, hit_counts.get(i)[0]);
                            locations[mate][genome][sequence].clear();
                        }
                    }
                }
            }
        }
        
        public void examine_location(int mate, int genome, int sequence, int ref_start){
            int start, stop, offset = 0, length, deletions;
            double identity = -1;
            String cigar;
            boolean forward = ref_start >= 0;
            single_hit curr_hit, similar_subject;
            if (!forward)
                ref_start = -ref_start - 1;
            reference.setLength(0);
            start = ref_start - ALIGNMENT_BOUND;
            stop = start + reads[mate].length() + 2 * ALIGNMENT_BOUND - 1;
            //System.out.println(ref_start + " " + sequence_length[genome][sequence]);
            if (start >= 0 && stop <= sequence_length[genome][sequence] - 1){
                number_of_hits.getAndIncrement();
                genome_scanner.get_sub_sequence(reference, genome, sequence, start, stop - start + 1, true);
                if ((similar_subject = find_similar_subject(reference, genome, mate)) != null){
                    cigar = similar_subject.cigar;
                    identity = similar_subject.identity;
                    offset = similar_subject.offset;
                    length = similar_subject.length;
                    deletions = similar_subject.deletions;
                    curr_hit = new single_hit(genome, sequence, identity, start, offset, length, deletions, forward, cigar, reference.toString());
                    hits[mate][genome].add(curr_hit);
                } else {
                    number_of_alignments.getAndIncrement();
                    bounded_aligner.align(forward?reads[mate].forward_seq:reads[mate].reverse_seq, reference); 
                    cigar = bounded_aligner.get_cigar();
                    identity = bounded_aligner.get_identity();
                    //System.out.println(cigar.toString());
                    offset = bounded_aligner.get_offset();
                    length = bounded_aligner.get_range_length();
                    deletions = bounded_aligner.get_deletions();
                    curr_hit = new single_hit(genome, sequence, identity, start, offset, length, deletions, forward, cigar, reference.toString());
                    //System.out.println(curr_hit.toString());
                    hits[mate][genome].add(curr_hit);
                    alignments[mate][genome].add(curr_hit);
                }
            } else if(ref_start >= 0 && ref_start + reads[mate].length() <= sequence_length[genome][sequence]){ // may map at the borders
                number_of_hits.getAndIncrement();
                genome_scanner.get_sub_sequence(reference, genome, sequence, ref_start, reads[mate].length(), true);
                if ((similar_subject = find_similar_subject(reference, genome, mate)) != null){
                    cigar = similar_subject.cigar;
                    identity = similar_subject.identity;
                    offset = similar_subject.offset;
                    length = similar_subject.length;
                    deletions = similar_subject.deletions;
                    curr_hit = new single_hit(genome, sequence, identity, ref_start, offset, length, deletions, forward, cigar, reference.toString());
                    hits[mate][genome].add(curr_hit);
                } else {
                    number_of_alignments.getAndIncrement();
                    aligner.align(forward?reads[mate].forward_seq:reads[mate].reverse_seq, reference); 
                    cigar = aligner.get_cigar();
                    identity = aligner.get_identity();
                    //System.out.println(cigar.toString());
                    offset = aligner.get_offset();
                    length = aligner.get_range_length();
                    deletions = aligner.get_deletions();
                    curr_hit = new single_hit(genome, sequence, identity, ref_start,  offset, length, deletions, forward, cigar, reference.toString());
                    //System.out.println(curr_hit.toString());
                    hits[mate][genome].add(curr_hit);
                    alignments[mate][genome].add(curr_hit);
                }
            }
        }

        public single_hit find_similar_subject(StringBuilder reference, int curr_genome, int mate) {
            single_hit similar_alignment = null;
            int genome;
            for (ListIterator<Integer> genome_itr = genome_numbers.listIterator();genome_itr.hasNext();){
                genome = genome_itr.next();
                if (genome > curr_genome)
                    break;
                Iterator<single_hit> itr = alignments[mate][genome].iterator();
                while (itr.hasNext()) {
                    similar_alignment = itr.next();
                    if (are_equal(similar_alignment.reference, reference))
                        return similar_alignment;
                }
            }
            return null;
        }  
        
        public boolean are_equal(String s1, StringBuilder s2){
            boolean are_equal = s1.length() == s2.length();
            for (int i = 0; are_equal && i < s1.length(); ++i)
                if (s1.charAt(i) != s2.charAt(i))
                    are_equal = false;
            return are_equal;
        }
        
        public void report_hits(){
            int genome;
            if (ALIGNMENT_MODE < 0){ // pan-genomic best
                for (ListIterator<Integer> genome_itr = genome_numbers.listIterator();genome_itr.hasNext();){
                    genome = genome_itr.next();
                    collect_hits(genome);
                }
                call_mode();
                clear_hits_list();
            } else if (ALIGNMENT_MODE > 0){ //genomic best
                for (ListIterator<Integer> genome_itr = genome_numbers.listIterator();genome_itr.hasNext();){
                    genome = genome_itr.next();
                    collect_hits(genome);
                    call_mode();
                    clear_hits_list();
                }
            }
            clear_alignments();
        }
        
        public void call_mode(){
            switch (Math.abs(ALIGNMENT_MODE)){
                case 1: // unique best hits
                    report_unique_hit();
                break;    
                case 2: // one best hits
                    report_one_hit();
                break;    
                case 3: // all best hits
                    report_all_hit(true);
                break;
                case 4: // all hits
                    report_all_hit(false);
            }            
        }
        
        public void collect_hits(int genome){
            if (!paired_end){
                single_hit h;
                boolean found = false;
                while (!hits[0][genome].isEmpty()){
                    h = hits[0][genome].remove();
                    if (h.identity <= MIN_IDENTITY || h.length < MIN_HIT_LENGTH || h.start + h.offset < 0 || h.start + h.offset + h.deletions + reads[0].length() >= sequence_length[h.genome][h.sequence])
                        continue;
                    single_hits.add(h);
                    found = true;
                }
                hits[0][genome].clear();
                if (!found)
                    single_hits.add(new single_hit(genome,0,-1,-1,0,0,0,true,null,null));
            } else {
                single_hit h1, h2, first_hit1, first_hit2;
                paired_hit h;
                Iterator<single_hit> itr;
                first_hit1 = hits[0][genome].peek();
                first_hit2 = hits[1][genome].peek();
                boolean found = false;
                while (!hits[0][genome].isEmpty()){
                    h1 = hits[0][genome].remove();
                    if (h1.identity <= MIN_IDENTITY || h1.length < MIN_HIT_LENGTH || h1.start + h1.offset < 0 || h1.start + h1.offset + h1.deletions + reads[0].length() >= sequence_length[h1.genome][h1.sequence])
                        continue;
                    for (itr = hits[1][genome].iterator(); itr.hasNext(); ){
                        h2 = itr.next();
                        if (h2.identity <= MIN_IDENTITY || h2.length < MIN_HIT_LENGTH || h2.start + h2.offset < 0 || h2.start + h2.offset + h2.deletions + reads[1].length() >= sequence_length[h2.genome][h2.sequence])
                            continue;
                        h = new paired_hit(fragment_length(h1, h2), h1.genome, h1.sequence, h1.identity, h1.start, h1.offset, h1.length, h1.deletions, h1.forward, h1.cigar, h1.reference,
                                           h2.genome, h2.sequence, h2.identity, h2.start, h2.offset, h2.length, h2.deletions, h2.forward, h2.cigar, h2.reference);
                        found = true;
                        paired_hits.add(h);
                    }
                }
                hits[0][genome].clear();// is already empty
                hits[1][genome].clear();
                if (!found){
                    h1 = first_hit1;
                    h2 = first_hit2;
                    if (hits[0][genome].isEmpty() || h1.identity <= MIN_IDENTITY){
                        if (hits[1][genome].isEmpty() || h2.identity <= MIN_IDENTITY)
                            paired_hits.add(new paired_hit(Integer.MAX_VALUE, genome,0,-1,-1,0,0,0,true,null,null,genome,0,-1,-1,0,0,0,true,null,null));
                        else
                            paired_hits.add(new paired_hit(Integer.MAX_VALUE, genome,0,-1,-1,0,0,0,true,null,null,h2.genome, h2.sequence, h2.identity, h2.start, h2.offset, h2.length, h2.deletions, h2.forward, h2.cigar, h2.reference));
                    } else {
                        if (hits[1][genome].isEmpty() || h2.identity <= MIN_IDENTITY)
                            paired_hits.add(new paired_hit(Integer.MAX_VALUE, h1.genome, h1.sequence, h1.identity, h1.start, h1.offset, h1.length, h1.deletions, h1.forward, h1.cigar, h1.reference,genome,0,-1,-1,0,0,0,true,null,null));
                        else
                            paired_hits.add(new paired_hit(fragment_length(h1, h2), h1.genome, h1.sequence, h1.identity, h1.start, h1.offset, h1.length, h1.deletions, h1.forward, h1.cigar, h1.reference,h2.genome, h2.sequence, h2.identity, h2.start, h2.offset, h2.length, h2.deletions, h2.forward, h2.cigar, h2.reference));
                    }
                }
            }
        }

        public void report_unique_hit(){
            if (!paired_end){
                single_hit h, best_hit;
                h = single_hits.remove(); 
                best_hit = new single_hit(h);
                if (h.start != -1){
                    if(!single_hits.isEmpty()){
                        h = single_hits.remove();
                        if (h.identity != best_hit.identity){
                            write_single_sam_record(best_hit, 0);
                            unique[best_hit.genome]++;
                        } else
                            shared[best_hit.genome]++;
                    } else {
                        write_single_sam_record(best_hit, 0);
                        unique[best_hit.genome]++;
                    }
                } else {
                    unmapped[best_hit.genome]++;
                    write_single_sam_record(best_hit, 4);
                }
            } else {
                paired_hit h, best_hit;
                h = paired_hits.remove(); 
                best_hit = new paired_hit(h);
                if (best_hit.get_max_start() != -1){
                    if(!paired_hits.isEmpty()){
                        h = paired_hits.remove();
                        if (h.get_identity() != best_hit.get_identity()){
                            if (best_hit.get_min_start() != -1){
                                write_paired_sam_record(best_hit, new int[]{0, 0});
                                unique[best_hit.genome1]++;
                                unique[best_hit.genome2]++;
                            } else if (best_hit.start1 != -1){
                                write_paired_sam_record(best_hit, new int[]{0, 4});
                                unique[best_hit.genome1]++;
                                unmapped[best_hit.genome2]++;
                            } else if (best_hit.start2 != -1){
                                write_paired_sam_record(best_hit, new int[]{4, 0});
                                unique[best_hit.genome2]++;
                                unmapped[best_hit.genome1]++;
                            }
                        } else {
                            if (best_hit.get_min_start() != -1){
                                write_paired_sam_record(best_hit, new int[]{0, 0});
                                shared[best_hit.genome1]++;
                                shared[best_hit.genome2]++;
                            } else if (best_hit.start1 != -1){
                                write_paired_sam_record(best_hit, new int[]{0, 4});
                                shared[best_hit.genome1]++;
                                unmapped[best_hit.genome2]++;
                            } else if (best_hit.start2 != -1){
                                write_paired_sam_record(best_hit, new int[]{4, 0});
                                shared[best_hit.genome2]++;
                                unmapped[best_hit.genome1]++;
                            }
                        }
                    } else {
                        if (best_hit.get_min_start() != -1){
                            write_paired_sam_record(best_hit, new int[]{0, 0});
                            unique[best_hit.genome1]++;
                            unique[best_hit.genome2]++;
                        } else if (best_hit.start1 != -1){
                            write_paired_sam_record(best_hit, new int[]{0, 4});
                            unique[best_hit.genome1]++;
                            unmapped[best_hit.genome2]++;
                        } else if (best_hit.start2 != -1){
                            write_paired_sam_record(best_hit, new int[]{4, 0});
                            unique[best_hit.genome2]++;
                            unmapped[best_hit.genome1]++;
                        }
                    }
                } else {
                    unmapped[best_hit.genome1]++;
                    unmapped[best_hit.genome2]++;
                    write_paired_sam_record(best_hit, new int[]{4, 4});
                }
            }
        }

        public void report_one_hit(){
            if (!paired_end){
                single_hit h, best_hit;
                int count;
                double rnd, freq = 0, sum_freq = 0;
                best_hit = single_hits.peek();
                for (count = 0; !single_hits.isEmpty(); ++count){
                    h = single_hits.remove();
                    if (h.start == -1 || h.identity != best_hit.identity)
                        break;
                    sum_freq += raw_abundance[h.genome] / genome_sizes[h.genome];
                    single_hits_2.add(h);
                }
                rnd = rand.nextDouble();
                while(!single_hits_2.isEmpty()){
                    best_hit = single_hits_2.remove();
                    freq += raw_abundance[best_hit.genome] / genome_sizes[best_hit.genome];
                    if (rnd < freq / sum_freq)
                        break;
                }
                if (best_hit.start != -1){
                    if (count == 1)
                        unique[best_hit.genome]++;
                    else
                        shared[best_hit.genome]++;
                    write_single_sam_record(best_hit, 0);
                } else {
                       unmapped[best_hit.genome]++;
                       write_single_sam_record(best_hit, 4);
                }   
            } else { // paired
                paired_hit h, best_hit;
                int count;
                double rnd, freq = 0, sum_freq = 0;
                best_hit = paired_hits.peek();
                for (count = 0; !paired_hits.isEmpty(); ++count){
                    h = paired_hits.remove();
                    if (h.get_max_start() == -1 || h.get_identity() != best_hit.get_identity())
                        break;
                    sum_freq += raw_abundance[h.genome1] / genome_sizes[h.genome1];
                    paired_hits_2.add(h);
                }
                rnd = rand.nextDouble();
                while(!paired_hits_2.isEmpty()){
                    best_hit = paired_hits_2.remove();
                    freq += raw_abundance[best_hit.genome1] / genome_sizes[best_hit.genome1];
                    if (rnd < freq / sum_freq)
                        break;
                }
                if (best_hit.get_max_start() != -1){
                    if (best_hit.get_min_start() != -1){
                        write_paired_sam_record(best_hit, new int[]{0, 0});
                        if (count == 1){
                            unique[best_hit.genome1]++;
                            unique[best_hit.genome2]++;
                        } else {
                            shared[best_hit.genome1]++;
                            shared[best_hit.genome2]++;
                        }
                    } else if (best_hit.start1 != -1){
                        write_paired_sam_record(best_hit, new int[]{0, 4});
                        unmapped[best_hit.genome2]++;
                        if (count == 1)
                            unique[best_hit.genome1]++;
                        else 
                            shared[best_hit.genome1]++;
                    } else if (best_hit.start2 != -1){
                        write_paired_sam_record(best_hit, new int[]{4, 0});
                        unmapped[best_hit.genome1]++;
                        if (count == 1)
                            unique[best_hit.genome2]++;
                        else 
                            shared[best_hit.genome2]++;
                    }
                } else { 
                       write_paired_sam_record(best_hit, new int[]{4, 4});
                        unmapped[best_hit.genome1]++;
                        unmapped[best_hit.genome2]++;
                }
            }
        }
        
        public void report_all_hit(boolean best){
            if (!paired_end){
                single_hit h, best_hit;
                int count, prev_genome;
                best_hit = single_hits.peek();
                for (count = 0; !single_hits.isEmpty(); ++count){
                    h = single_hits.remove();
                    if (h.start == -1 || (best && h.identity != best_hit.identity))
                        break;
                    single_hits_2.add(h); 
                }
                if (best_hit.start != -1){
                    prev_genome = -1;
                    while (!single_hits_2.isEmpty()){
                        best_hit = single_hits_2.remove();
                        if (count == 1)
                            unique[best_hit.genome]++;
                        else
                            shared[best_hit.genome]++;
                        write_single_sam_record(best_hit, best_hit.genome == prev_genome ? 256 : 0);
                        prev_genome = best_hit.genome;
                    }
                } else {
                       unmapped[best_hit.genome]++;
                       write_single_sam_record(best_hit, 4);
                } 
            } else {
                paired_hit h, best_hit;
                int count, prev_genome;
                best_hit = paired_hits.peek();
                for (count = 0; !paired_hits.isEmpty(); ++count){
                    h = paired_hits.remove();
                    if (h.get_max_start() == -1 || (best && h.get_identity() != best_hit.get_identity()))
                        break;
                    paired_hits_2.add(h);
                }
                if (best_hit.get_max_start() != -1){
                    prev_genome = -1;
                    while (!paired_hits_2.isEmpty()){
                        best_hit = paired_hits_2.remove();
                        if (best_hit.get_min_start() != -1){
                            write_paired_sam_record(best_hit, new int[]{best_hit.genome1 == prev_genome ? 256 : 0, best_hit.genome2 == prev_genome ? 256 : 0});
                            if (count == 1){
                                unique[best_hit.genome1]++;
                                unique[best_hit.genome2]++;
                            } else {
                                shared[best_hit.genome1]++;
                                shared[best_hit.genome2]++;
                            }
                        } else if (best_hit.start1 != -1){
                            write_paired_sam_record(best_hit, new int[]{best_hit.genome1 == prev_genome ? 256 : 0, 4});
                            unmapped[best_hit.genome2]++;
                            if (count == 1)
                                unique[best_hit.genome1]++;
                            else 
                                shared[best_hit.genome1]++;
                        } else if (best_hit.start2 != -1){
                            write_paired_sam_record(best_hit, new int[]{4, best_hit.genome2 == prev_genome ? 256 : 0});
                            unmapped[best_hit.genome1]++;
                            if (count == 1)
                                unique[best_hit.genome2]++;
                            else 
                                shared[best_hit.genome2]++;
                        }
                        prev_genome = best_hit.genome1;
                    }
                } else { 
                       write_paired_sam_record(best_hit, new int[]{4, 4});
                        unmapped[best_hit.genome1]++;
                        unmapped[best_hit.genome2]++;
                }
                
            }
        }
        

        public void clear_hits_list(){
            if (!paired_end) {
                single_hits.clear();
                single_hits_2.clear();
            } else {
                paired_hits.clear();
                paired_hits_2.clear();
            }
        }
        
        public void clear_alignments(){    
            for (ListIterator<Integer> genome_itr = genome_numbers.listIterator();genome_itr.hasNext();)
                alignments[0][genome_itr.next()].clear();
            if (paired_end)
                for (ListIterator<Integer> genome_itr = genome_numbers.listIterator();genome_itr.hasNext();)
                    alignments[1][genome_itr.next()].clear();
        }
        

        public void write_single_sam_record(single_hit h, int flag){
            SAMRecord sam_record = new SAMRecord(null);
            if (h.start == -1){
                sam_record.setReadName(reads[0].name.toString());
                sam_record.setFlags(flag);
                sam_record.setReadString(reads[0].forward_seq.toString());
            } else {
                flag |= h.forward?0:16; // direction
                sam_record.setReferenceName(sequence_titles[h.genome][h.sequence].split("\\s")[0]);
                sam_record.setFlags(flag);
                sam_record.setReadName(reads[0].name.toString());
                sam_record.setAlignmentStart(h.start + h.offset + 1);
                sam_record.setMappingQuality((int)Math.round(-10 * Math.log10(1 - (h.identity - 1) / 100.0)) + 1);
                sam_record.setCigarString(h.cigar);
                sam_record.setReadString((h.forward?reads[0].forward_seq:reads[0].reverse_seq).toString());
            }
            synchronized(sam_headers[h.genome]){
                synchronized(sam_writers[h.genome]){
                    sam_record.setHeader(sam_headers[h.genome]);
                    sam_writers[h.genome].addAlignment(sam_record);
                }
            }
        }
        
        public void write_paired_sam_record(paired_hit h, int[] flag){
        // the first segment    
            SAMRecord sam_record = new SAMRecord(null);
            int position1 = h.start1 + h.offset1 + 1;
            int position2 = h.start2 + h.offset2 + 1;;
            flag[0] |= 64; 
            if (h.start1 == -1){
                sam_record.setReadName(reads[0].name.toString());
                sam_record.setFlags(flag[0]);
                sam_record.setReadString(reads[0].forward_seq.toString());
            } else {
                flag[0] |= h.forward1?0:16; // SEQ being reverse complemented
                flag[0] |= !h.forward2?32:0; // SEQ being reverse complemented
            //qname
                sam_record.setReadName(reads[0].name.toString()); 
            //flag    
                sam_record.setFlags(flag[0]); 
            //rname    
                sam_record.setReferenceName(sequence_titles[h.genome1][h.sequence1].split("\\s")[0]);
            //pos    
                sam_record.setAlignmentStart(position1);
            //mapq    
                sam_record.setMappingQuality((int)Math.round(-10 * Math.log10(1 - (h.identity1/100.0 - 0.01))) + 1);
                //sam_record.setMappingQuality((int)Math.round(h.identity1));
            //cigar    
                sam_record.setCigarString(h.cigar1);
            //rnext
                if (h.start2 == -1)
                    sam_record.setMateReferenceName("*");
                else 
                    sam_record.setMateReferenceName(sequence_titles[h.genome2][h.sequence2].split("\\s")[0]);
            //pnext 
                sam_record.setMateAlignmentStart(position2);
            //tlen
                if (h.start2 == -1)
                    sam_record.setInferredInsertSize(0);
                else if (position2 > position1)
                    sam_record.setInferredInsertSize(position2 - position1 + reads[0].forward_seq.length());
                else
                    sam_record.setInferredInsertSize(-(position1 - position2 + reads[1].forward_seq.length()));
            //seq    
                sam_record.setReadString((h.forward1?reads[0].forward_seq:reads[0].reverse_seq).toString());
            //qual    
            }
            synchronized(sam_headers[h.genome1]){
                synchronized(sam_writers[h.genome1]){
                    sam_record.setHeader(sam_headers[h.genome1]);
                    sam_writers[h.genome1].addAlignment(sam_record);
                }
            }
        // the second segment    
            sam_record = new SAMRecord(null);
            flag[1] |= 128; 
            if (h.start2 == -1){ // not mapped
                sam_record.setReadName(reads[1].name.toString());
                sam_record.setFlags(flag[1]);
                sam_record.setReadString(reads[1].forward_seq.toString());
            } else {
                flag[1] |= h.forward2?0:16; // SEQ being reverse complemented
                flag[1] |= !h.forward1?32:0; // SEQ being reverse complemented
            //qname
                sam_record.setReadName(reads[1].name.toString()); 
            //flag    
                sam_record.setFlags(flag[1]); 
            //rname    
                sam_record.setReferenceName(sequence_titles[h.genome2][h.sequence2].split("\\s")[0]);
            //pos    
                sam_record.setAlignmentStart(position2);
            //mapq    
                sam_record.setMappingQuality((int)Math.round(-10 * Math.log10(1 - (h.identity2/100.0 - 0.01))) + 1);
                //sam_record.setMappingQuality((int)Math.round(h.identity2));
            //cigar    
                sam_record.setCigarString(h.cigar2);
            //rnext
                if (h.start1 == -1)
                    sam_record.setMateReferenceName("*");
                else 
                    sam_record.setMateReferenceName(sequence_titles[h.genome1][h.sequence1].split("\\s")[0]);
            //pnext    
                sam_record.setMateAlignmentStart(position1);
            //tlen
                if (h.start1 == -1)
                    sam_record.setInferredInsertSize(0);
                else if (position2 > position1)
                    sam_record.setInferredInsertSize(-(position2 - position1 + reads[0].forward_seq.length()));
                else
                    sam_record.setInferredInsertSize(position1 - position2 + reads[1].forward_seq.length());
            //seq    
                sam_record.setReadString((h.forward2?reads[1].forward_seq:reads[1].reverse_seq).toString());
            //qual    
            }
            synchronized(sam_headers[h.genome2]){
                synchronized(sam_writers[h.genome2]){
                    sam_record.setHeader(sam_headers[h.genome2]);
                    sam_writers[h.genome2].addAlignment(sam_record);
                }
            }
        }
        

        public int fragment_length(single_hit h1, single_hit h2){
            int frag_len;
            int position1, position2;
            if (h1.genome == h2.genome && h1.sequence == h2.sequence && h1.forward == !h2.forward){
                position1 = h1.start + h1.offset;
                position2 = h2.start + h2.offset;
                if (h1.forward)
                    frag_len = position2 + reads[1].forward_seq.length() - position1;
                else
                    frag_len = position1 + reads[0].forward_seq.length() - position2;
                if (frag_len >= Math.max(reads[0].length(), reads[1].length()) && frag_len < MAX_FRAGMENT_LENGTH)
                    return frag_len;
            }
            return Integer.MAX_VALUE;
        }
    }   
    
    public void map_reads() {
        int i, number, genome, n = 0;
        Node pangenome_node;
        LinkedList<Integer> genome_numbers, gn;
        String line;
        BufferedReader in;
        BufferedWriter out;
        int total_unique, total_mapped, total_unmapped;
        boolean paired;
        File theDir;
        Scanner s;
        String str;
        s = new Scanner(System.in);
        if (PATH_TO_THE_FIRST_SRA == null){
            System.out.println("PATH_TO_THE_FIRST_SRA is empty.");
            System.exit(1);
        }
        paired = PATH_TO_THE_SECOND_SRA != null;
        if (PATH_TO_THE_GENOME_NUMBERS_FILE == null){
            System.out.println("PATH_TO_THE_GENOME_NUMBERS_FILE is empty.");
            System.exit(1);
        }
        if (! new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH).exists()) {
            System.out.println("No database found in " + PATH_TO_THE_PANGENOME_DATABASE);
            System.exit(1);
        }
        if (OUTPUT_PATH.equals(""))
            OUTPUT_PATH = PATH_TO_THE_PANGENOME_DATABASE;
        if (! new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH).exists()) {
            System.out.println("No graph database found in " + PATH_TO_THE_PANGENOME_DATABASE);
            System.exit(1);
        }
        if (! new File(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH).exists()) {
            System.out.println("No index database found in " + PATH_TO_THE_PANGENOME_DATABASE);
            System.exit(1);
        }
        if (! new File(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH).exists()) {
            System.out.println("No genome database found in " + PATH_TO_THE_PANGENOME_DATABASE);
            System.out.println("Do you want to reconstruct it from the graph database [y/n]? ");
            str = s.nextLine().toLowerCase();
            while (!str.equals("y") && !str.equals("n")){
            System.out.println("Do you want to reconstruct it from the graph database [y/n]? ");
                 str = s.nextLine().toLowerCase();
            }
            if (str.equals("y")){
                rebuild_genome_database();
            } else {
                System.out.println("Exiting the program...");
                System.exit(1);  
            }
        } else
            genomeDb = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH);
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
        registerShutdownHook(graphDb);
        startTime = System.currentTimeMillis();
        try (Transaction tx = graphDb.beginTx()) {
            pangenome_node = graphDb.findNodes(pangenome_label).next();
            if (pangenome_node == null) {
                System.out.println("Can not locate database node!");
                System.exit(1);
            }
            tx.success();
        }
        num_shared_mapping = new AtomicInteger[genomeDb.num_genomes + 1];
        num_unique_mapping = new AtomicInteger[genomeDb.num_genomes + 1];
        num_unmapped = new AtomicInteger[genomeDb.num_genomes + 1];
        read_number = new AtomicInteger(1);
        number_of_alignments = new AtomicInteger(0);
        number_of_hits = new AtomicInteger(0);
        genome_numbers = new LinkedList();
        try {
            in = new BufferedReader(new FileReader(PATH_TO_THE_GENOME_NUMBERS_FILE));
            for (n = 0; in.ready(); ){
                line = in.readLine().trim();
                if (!line.equals("")){
                    number = Integer.parseInt(line);
                    if (number > 0 && number <= genomeDb.num_genomes){
                        genome_numbers.add(number);
                        num_shared_mapping[number] = new AtomicInteger(0);
                        num_unique_mapping[number] = new AtomicInteger(0);
                        num_unmapped[number] = new AtomicInteger(0);
                        ++n;
                    } else {
                        System.err.println("Genome " + number + "is not in the database");
                    }
                }
            }
            in.close();
        } catch (Exception ex){
            System.err.println("Error in reading genome numbers");
        }
        SAMFileWriter[] sams = new SAMFileWriter[genomeDb.num_genomes + 1];
        SAMFileHeader[] headers = new SAMFileHeader[genomeDb.num_genomes + 1];
        for (ListIterator<Integer> genome_itr = genome_numbers.listIterator();genome_itr.hasNext();){
            genome = genome_itr.next();
            headers[genome] = new SAMFileHeader();
            for (i = 1; i <= genomeDb.num_sequences[genome]; ++i)
                headers[genome].addSequence(new SAMSequenceRecord(genomeDb.sequence_titles[genome][i].split("\\s")[0], (int)genomeDb.sequence_length[genome][i]));
            headers[genome].addProgramRecord(new SAMProgramRecord("PanTools"));
            if (BAMFORMAT)
                sams[genome] = new SAMFileWriterFactory().makeBAMWriter(headers[genome], false, 
                        new File(OUTPUT_PATH + "/pantools_" + genome + ".bam"));
            else
                sams[genome] = new SAMFileWriterFactory().makeSAMWriter(headers[genome], false, 
                        new File(OUTPUT_PATH + "/pantools_" + genome + ".sam"));
        }
        SequenceDatabase[] sequencingDb = new SequenceDatabase[2];
        theDir = new File(OUTPUT_PATH + "/sra1");
        if (theDir.exists()) {
            System.out.println("A sequencing databases already exists at " + OUTPUT_PATH + "/sra1" + ".");
            System.out.println("Do you want to reuse it [y/n]? ");
            str = s.nextLine().toLowerCase();
            while (!str.equals("y") && !str.equals("n")){
                 System.out.println("Do you want to reuse it [y/n]? ");
                 str = s.nextLine().toLowerCase();
            }
            if (str.equals("y"))
               sequencingDb[0] = new SequenceDatabase(OUTPUT_PATH + "/sra1");
            else{
                try {
                    FileUtils.deleteRecursively(new File(OUTPUT_PATH + "/sra1"));
                } catch (IOException ioe) {
                    System.out.println("Failed to delete the sequencing database");
                    System.exit(1);  
                }
                sequencingDb[0] = new SequenceDatabase(OUTPUT_PATH + "/sra1", PATH_TO_THE_FIRST_SRA);
            } 
        } else
            sequencingDb[0] = new SequenceDatabase(OUTPUT_PATH + "/sra1", PATH_TO_THE_FIRST_SRA);
        if (paired){
            theDir = new File(OUTPUT_PATH + "/sra2");
            if (theDir.exists()) {
                System.out.println("A sequencing databases already exists at " + OUTPUT_PATH + "/sra2" + ".");
                System.out.println("Do you want to reuse it [y/n]? ");
                str = s.nextLine().toLowerCase();
                while (!str.equals("y") && !str.equals("n")){
                     System.out.println("Do you want to reuse it [y/n]? ");
                     str = s.nextLine().toLowerCase();
                }
                if (str.equals("y"))
                   sequencingDb[1] = new SequenceDatabase(OUTPUT_PATH + "/sra2");
                else{
                    try {
                        FileUtils.deleteRecursively(new File(OUTPUT_PATH + "/sra2"));
                    } catch (IOException ioe) {
                        System.out.println("Failed to delete the sequencing database");
                        System.exit(1);  
                    }
                    sequencingDb[1] = new SequenceDatabase(OUTPUT_PATH + "/sra2", PATH_TO_THE_SECOND_SRA);
                } 
            } else
                sequencingDb[1] = new SequenceDatabase(OUTPUT_PATH + "/sra2", PATH_TO_THE_SECOND_SRA);
            if (sequencingDb[0].num_sequences[1] != sequencingDb[1].num_sequences[1]){
                System.err.println("Paired-end libraries contain different number of reads.");
                new File(OUTPUT_PATH + "/sra1/sequences.db").delete();
                new File(OUTPUT_PATH + "/sra2/sequences.db").delete();
                new File(OUTPUT_PATH + "/sra1/sequences.info").delete();
                new File(OUTPUT_PATH + "/sra2/sequences.info").delete();
                new File(OUTPUT_PATH + "/sra1").delete();
                new File(OUTPUT_PATH + "/sra2").delete();
                System.exit(1);
            }
        }

        System.out.println("Mapping " + sequencingDb[0].num_sequences[1] + (paired?" paired-end ": " ") + "reads on " + n + " genome(s) :");
        if (n > 0){
            System.out.print("0..................................................100\n ");
            try{
                ExecutorService es = Executors.newFixedThreadPool(THREADS);
                for (i = 0; i < THREADS; ++i){
                    gn = new LinkedList();
                    for (ListIterator<Integer> genome_itr = genome_numbers.listIterator();genome_itr.hasNext();)
                         gn.add((int)genome_itr.next());
                    es.execute(new Map(i, gn, paired, sams, headers, sequencingDb));
                }
                es.shutdown();
                es.awaitTermination(10, TimeUnit.DAYS);        
            } catch (InterruptedException e){

            }
            
            System.out.println("\nNumber_of_reads = " + sequencingDb[0].num_sequences[1] * (paired ? 2 : 1));
            System.out.println("Number_of_hits = " + number_of_hits);
            System.out.println("Number_of_alignments = " + number_of_alignments);

            total_unique = total_mapped = total_unmapped = 0;
            try{
                out = new BufferedWriter(new FileWriter(OUTPUT_PATH + "/mapping_summary.txt"));
                out.write("\nGenome\tTotal\tUnique\tUnmapped\n");
                System.out.print("\nGenome\tTotal\tUnique\tUnmapped\n");
                for (ListIterator<Integer> genome_itr = genome_numbers.listIterator();genome_itr.hasNext();){
                    genome = genome_itr.next();
                    total_unique += num_unique_mapping[genome].get();
                    total_mapped += num_shared_mapping[genome].get() +
                                    num_unique_mapping[genome].get();
                    total_unmapped += num_unmapped[genome].get();
                    out.write(genome + "\t" + (num_shared_mapping[genome].get() + num_unique_mapping[genome].get())
                                              + "\t" + num_unique_mapping[genome].get() 
                                              + "\t" + num_unmapped[genome].get() + "\n");
                    System.out.print(genome + "\t" + (num_shared_mapping[genome].get() + num_unique_mapping[genome].get())
                                              + "\t" + num_unique_mapping[genome].get() 
                                              + "\t" + num_unmapped[genome].get() + "\n");
                    sams[genome].close();
                }
                out.close();
            } catch (IOException ex){
                System.err.println(ex.getMessage());
            }
            System.out.println("\nTotal mapped:\t" + total_mapped + "\nUniquely mapped:\t" + total_unique + "\nUnmapped:\t" + total_unmapped);
        }
        genomeDb.close();
        graphDb.shutdown();
    }

    /**
     * Constructs a pangenome database from given genomes.
     * 
     * @param genome_paths_file Path to the FASTA genome files. 
     * @param pangenome_path Path to the database folder
     */  
    public void initialize_pangenome() {
        Node pangenome_node;
        create_database();
        try (Transaction tx = graphDb.beginTx()) {
            pangenome_node = graphDb.createNode(pangenome_label);
            pangenome_node.setProperty("k_mer_size", K_SIZE);
            pangenome_node.setProperty("date", new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(new Date()));
            tx.success();
        }

        num_nodes = 0;
        num_edges = 0;
        num_bases = 0;
        num_degenerates = 0;
        construct_pangenome(pangenome_node);
        System.out.println("Number of kmers:   " + indexDb.length());
        System.out.println("Number of nodes:   " + num_nodes);
        System.out.println("Number of edges:   " + num_edges);
        System.out.println("Number of bases:   " + num_bases);
        System.out.println("Number of degenerate nodes:   " + num_degenerates);
        try (Transaction tx = graphDb.beginTx()) {
            pangenome_node.setProperty("k_mer_size", K_SIZE);
            pangenome_node.setProperty("num_k_mers", indexDb.length());
            pangenome_node.setProperty("num_nodes", num_nodes);
            pangenome_node.setProperty("num_degenerate_nodes", num_degenerates);
            pangenome_node.setProperty("num_edges", num_edges);
            pangenome_node.setProperty("num_genomes", genomeDb.num_genomes);
            pangenome_node.setProperty("num_bases", num_bases);
            tx.success();
        }
        graphDb.shutdown();
        genomeDb.close();
        indexDb.close();
        File directory = new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH);
        for (File f : directory.listFiles()) {
            if (f.getName().startsWith("neostore.transaction.db.")) {
                f.delete();
            }
        }
        System.out.println("graph.db size: " + getFolderSize(new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH)) + " MB");
        System.out.println("index.db size: " + getFolderSize(new File(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH)) + " MB");
        System.out.println("genome.db size: " + getFolderSize(new File(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH)) + " MB");
    }
    
    /**
     * Adds new genomes to an available pangenome.
     * 
     * @param genome_paths_file Path to the FASTA genome files. 
     * @param pangenome_path Path to the database folder
     */
    public void add_genomes() {
        int previous_num_genomes;
        Node pangenome_node;
        Scanner s;
        String str;
        s = new Scanner(System.in);
        if (PATH_TO_THE_GENOMES_FILE == null){
            System.out.println("PATH_TO_THE_GENOMES_FILE is empty.");
            System.exit(1);
        }
        if (! new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH).exists()) {
            System.out.println("No graph database found in " + PATH_TO_THE_PANGENOME_DATABASE);
            System.exit(1);
        }
        if (! new File(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH).exists()) {
            System.out.println("No index database found in " + PATH_TO_THE_PANGENOME_DATABASE);
            System.exit(1);
        }
        if (! new File(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH).exists()) {
            System.out.println("No genome database found in " + PATH_TO_THE_PANGENOME_DATABASE);
            System.out.println("Do you want to reconstruct it from the graph database [y/n]? ");
            str = s.nextLine().toLowerCase();
            while (!str.equals("y") && !str.equals("n")){
            System.out.println("Do you want to reconstruct it from the graph database [y/n]? ");
                 str = s.nextLine().toLowerCase();
            }
            if (str.equals("y")){
                rebuild_genome_database();
            } else {
                System.out.println("Exiting the program...");
                System.exit(1);  
            }
        } else
            genomeDb = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH);
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
        registerShutdownHook(graphDb);
        startTime = System.currentTimeMillis();
        try (Transaction tx = graphDb.beginTx()) {
            pangenome_node = graphDb.findNodes(pangenome_label).next();
            if (pangenome_node == null) {
                System.out.println("Can not locate database node!");
                System.exit(1);
            }
        // Reads the properties of the pangenome    
            K_SIZE = (int) pangenome_node.getProperty("k_mer_size");
            num_nodes = (long) pangenome_node.getProperty("num_nodes");
            num_edges = (long) pangenome_node.getProperty("num_edges");
            num_degenerates = (int) pangenome_node.getProperty("num_degenerate_nodes");
            num_bases = 0;
            previous_num_genomes = (int) pangenome_node.getProperty("num_genomes");
            tx.success();
        }
        genomeDb.add_genomes(PATH_TO_THE_GENOMES_FILE);
        indexDb = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH, PATH_TO_THE_GENOMES_FILE, genomeDb, graphDb, previous_num_genomes);
        genome_scanner = new SequenceScanner(genomeDb, previous_num_genomes + 1, 1, K_SIZE, indexDb.get_pre_len());
    // the sequences should be dropped out as they will change and add_sequence_properties() function will rebuild them.    
        drop_nodes_property("sequence");
    // the edge colors should be dropped out as they will change and localize_nodes() function will rebuild them again.    
        drop_edges_colors();
        construct_pangenome(pangenome_node);
        System.out.println("Number of kmers:   " + indexDb.length());
        System.out.println("Number of nodes:   " + num_nodes);
        System.out.println("Number of edges:   " + num_edges);
        System.out.println("Number of bases:   " + num_bases);
        System.out.println("Number of degenerate nodes:   " + num_degenerates);
        try (Transaction tx = graphDb.beginTx()) {
            pangenome_node.setProperty("k_mer_size", K_SIZE);
            pangenome_node.setProperty("num_k_mers", indexDb.length());
            pangenome_node.setProperty("num_nodes", num_nodes);
            pangenome_node.setProperty("num_degenerate_nodes", num_degenerates);
            pangenome_node.setProperty("num_edges", num_edges);
            pangenome_node.setProperty("num_genomes", genomeDb.num_genomes);
            pangenome_node.setProperty("num_bases", num_bases);
            tx.success();
        }
        graphDb.shutdown();
        genomeDb.close();
        indexDb.close();
        File directory = new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH);
        for (File f : directory.listFiles()) {
            if (f.getName().startsWith("neostore.transaction.db.")) {
                f.delete();
            }
        }
        System.out.println("graph.db size: " + getFolderSize(new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH)) + " MB");
        System.out.println("index.db size: " + getFolderSize(new File(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH)) + " MB");
        System.out.println("genome.db size: " + getFolderSize(new File(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH)) + " MB");
    }

    public SequenceDatabase rebuild_genome_database(){
    // read genomes information from the graph and rebuild the genomes database
        int j, len;
        long byte_number = 0;
        int[] address = new int[4];
        Node start, seq_node;
        SequenceDatabase genomeDb = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH, graphDb);
        StringBuilder seq = new StringBuilder();
        for (address[0] = 1; address[0] <= genomeDb.num_genomes; ++address[0]) {
            for (address[1] = 1; address[1] <= genomeDb.num_sequences[address[0]]; ++address[1]) {
                seq_node = graphDb.findNode(sequence_label, "identifier", address[0] + "_" + address[1]);
                start = seq_node.getRelationships(Direction.OUTGOING).iterator().next().getEndNode();
                address[2] = 1;
                address[3] = (int) genomeDb.sequence_length[address[0]][address[1]];
                extract_sequence(seq, new IndexPointer(start.getId(), true, 0, -1l), address);
                len = seq.length();
                if (len % 2 == 1) {
                    --len;
                }
                for (j = 0; j < len; j += 2, ++byte_number) {
                    genomeDb.genomes_buff[(int) (byte_number / genomeDb.MAX_BYTE_COUNT)].put((byte) ((genomeDb.binary[seq.charAt(j)] << 4) | genomeDb.binary[seq.charAt(j + 1)]));
                }
                if (len == seq.length() - 1) {
                    genomeDb.genomes_buff[(int) (byte_number / genomeDb.MAX_BYTE_COUNT)].put((byte) (genomeDb.binary[seq.charAt(len)] << 4));
                    ++byte_number;
                }
            }
        }
        return genomeDb;
    }
    
    public void remove_genomes() {
        
    }
        
    public void create_database(){
        File theDir;
        Scanner s;
        String str;
        s = new Scanner(System.in);
        startTime = System.currentTimeMillis();
        if (PATH_TO_THE_PANGENOME_DATABASE == null){
            System.out.println("PATH_TO_THE_PANGENOME_DATABASE is empty.");
            System.exit(1);
        }
        theDir = new File(PATH_TO_THE_PANGENOME_DATABASE);
        if (theDir.exists()) {
        // creating a graph database
            theDir = new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH);
            if (theDir.exists()) {
                System.out.println("A graph databases already exists at " + PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH + ".");
                System.out.println("Do you want to remove it [y/n]? ");
                str = s.nextLine().toLowerCase();
                while (!str.equals("y") && !str.equals("n")){
                     System.out.println("Do you want to remove it [y/n]? ");
                     str = s.nextLine().toLowerCase();
                }
                if (str.equals("y")){
                    try {
                        FileUtils.deleteRecursively(new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH));
                    } catch (IOException ioe) {
                        System.out.println("Failed to delete the graph database");
                        System.exit(1);  
                    }
                } else {
                    System.out.println("Exiting the program...");
                    System.exit(1);  
                }
            }
            graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();  
            registerShutdownHook(graphDb);
        // connecting to genome database
            theDir = new File(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH);
            if (theDir.exists()) {
                System.out.println("A genome databases already exists at " + PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH + ".");
                System.out.println("Do you want to reuse it [y/n]? ");
                str = s.nextLine().toLowerCase();
                while (!str.equals("y") && !str.equals("n")){
                     System.out.println("Do you want to reuse it [y/n]? ");
                     str = s.nextLine().toLowerCase();
                }
                if (str.equals("y"))
                    genomeDb = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH);
                else {
                    try {
                        FileUtils.deleteRecursively(new File(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH));
                    } catch (IOException ioe) {
                        System.out.println("Failed to delete the genome database!");
                        System.exit(1);  
                    }                    
                    genomeDb = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH, PATH_TO_THE_GENOMES_FILE);
                }
            } else
                genomeDb = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH, PATH_TO_THE_GENOMES_FILE);
        // connecting to index database
            theDir = new File(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH);
            if (theDir.exists()) {
                System.out.println("An index databases already exists at " + PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH + ".");
                System.out.println("Do you want to reuse it [y/n]? ");
                str = s.nextLine().toLowerCase();
                while (!str.equals("y") && !str.equals("n")){
                     System.out.println("Do you want to reuse it [y/n]? ");
                     str = s.nextLine().toLowerCase();
                }
                if (str.equals("y")){
                    indexDb = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH, "sorted");
                    indexDb.clear_pointer_database();
                } else {
                    try {
                        FileUtils.deleteRecursively(new File(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH));
                    } catch (IOException ioe) {
                        System.out.println("Failed to delete the index database!");
                        System.exit(1);  
                    }                    
                }
                indexDb = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH, PATH_TO_THE_GENOMES_FILE, genomeDb, K_SIZE);
            } else
                indexDb = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH, PATH_TO_THE_GENOMES_FILE, genomeDb, K_SIZE);
            K_SIZE = indexDb.get_K();
            genome_scanner = new SequenceScanner(genomeDb, 1, 1, K_SIZE, indexDb.get_pre_len());
        } else {
            try {
                if (PATH_TO_THE_GENOMES_FILE == null){
                    System.out.println("PATH_TO_THE_GENOMES_FILE is empty.");
                    System.exit(1);
                }                
                System.out.println("Creating a database at " + PATH_TO_THE_PANGENOME_DATABASE);
                theDir.mkdir();
                graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH))
                        .setConfig(keep_logical_logs, "4 files").newGraphDatabase();  
                registerShutdownHook(graphDb);
                genomeDb = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH, PATH_TO_THE_GENOMES_FILE);
                indexDb = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH, PATH_TO_THE_GENOMES_FILE, genomeDb, K_SIZE);
                K_SIZE = indexDb.get_K();
                genome_scanner = new SequenceScanner(genomeDb, 1, 1, K_SIZE, indexDb.get_pre_len());
         } catch (SecurityException se) {
                System.out.println("Failed to create directory " + PATH_TO_THE_PANGENOME_DATABASE);
                System.exit(1);
            }
        }
        
    }
        
    /**
     * Retrieves the sequence of a number of genomic regions from the pangenome and stores them in a FASTA file.
     * 
     * @param region_records_file a text file with lines containing genome number, sequence number, start and stop positions
     *        of the genomic regions seperated by one space.
     * @param pangenome_path Path to the database folder
     */
    public void retrieve_regions() {
        if (PATH_TO_THE_REGIONS_FILE == null){
            System.out.println("PATH_TO_THE_REGIONS_FILE is empty.");
            System.exit(1);
        }
        if (! new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH).exists()) {
            System.out.println("No database found in " + PATH_TO_THE_PANGENOME_DATABASE);
            System.exit(1);
        }
        String[] fields;
        String line, out_file_name;
        IndexPointer start_ptr;
        StringBuilder seq;
        int c, num_regions = 0, proper_regions = 0;
        int[] address = new int[4];
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
        registerShutdownHook(graphDb);
        try (Transaction tx = graphDb.beginTx()) {
            K_SIZE = (int) graphDb.findNodes(pangenome_label).next().getProperty("k_mer_size");
            tx.success();
        }
        seq = new StringBuilder();
        try {
            BufferedReader in = new BufferedReader(new FileReader(PATH_TO_THE_REGIONS_FILE));
            while (in.ready()) {
                line = in.readLine().trim();
                if (line.equals("")) {
                    continue;
                }
                ++num_regions;
            }
            in.close();
        } catch (IOException e) {
            System.out.println(e.getMessage());
            System.exit(1);
        }
        startTime = System.currentTimeMillis();
        genomeDb = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH);
        indexDb = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH, "sorted");
        genome_scanner = new SequenceScanner(genomeDb, 1, 1, K_SIZE, indexDb.get_pre_len());
        try (Transaction tx = graphDb.beginTx()) {
            try (BufferedReader in = new BufferedReader(new FileReader(PATH_TO_THE_REGIONS_FILE))) {
                fields = PATH_TO_THE_REGIONS_FILE.split("\\/");
                out_file_name = PATH_TO_THE_PANGENOME_DATABASE + "/" + fields[fields.length - 1] + ".fasta";
                BufferedWriter out = new BufferedWriter(new FileWriter(out_file_name));
                for (c = 0; in.ready();) {
                    line = in.readLine().trim();
                    if (line.equals("")) {
                        continue;
                    }
                    fields = line.trim().split("\\s");
                    address[0] = Integer.parseInt(fields[0]);
                    address[1] = Integer.parseInt(fields[1]);
                    address[2] = Integer.parseInt(fields[2]);
                    address[3] = Integer.parseInt(fields[3]);
                    if (address[0] <= genomeDb.num_genomes && address[1] <= genomeDb.num_sequences[address[0]] && address[2] >= 1 && address[3] <= genomeDb.sequence_length[address[0]][address[1]]){
                        start_ptr = locate(graphDb, address, K_SIZE);
                        proper_regions++;
                        //extract_sequence(seq, start_ptr, address, K);
                        out.write(">genome:" + address[0] + " sequence:" + address[1] + " from:" + address[2] + " to:" + address[3] + " length:" + (address[3] - address[2] + 1) + "\n");
                        address[2] -= 1;
                        address[3] -= 1;
                        seq.setLength(0);
                        genome_scanner.get_sub_sequence(seq, address, true);
                        write_fasta(out, seq.toString(), 70);
                        ++c;
                        //if (c % (num_regions / 100 + 1) == 0) 
                        //    System.out.print((long) c * 100 / num_regions + 1 + "%\r");
                    } else
                        System.out.println(line + "is not a proper coordinate!");
                }
                in.close();
                out.close();
                System.out.println(proper_regions + " out of " + num_regions + " genomic regions found and retrieved successfully (See " + out_file_name + ")");
            } catch (IOException ioe) {
                System.out.println("Failed to read file names!");
                System.exit(1);
            }
            tx.success();
        }
        graphDb.shutdown();
    }
    
    /**
     * Reconstructs all or some of the genomes in separated FASTA files.
     * 
     * @param genome_numbers_file A text file containing the genome numbers to be . 
     * @param pangenome_path Path to the database folder
     */
    public void retrieve_genomes() {
        if (PATH_TO_THE_GENOME_NUMBERS_FILE == null){
            System.out.println("PATH_TO_THE_GENOME_NUMBERS_FILE is empty.");
            System.exit(1);
        }
        if (! new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH).exists()) {
            System.out.println("No database found in " + PATH_TO_THE_PANGENOME_DATABASE);
            System.exit(1);
        }
        BufferedReader in;
        BufferedWriter out;
        IndexPointer start;
        String genome_number;
        int[] address;
        StringBuilder seq;
        Node pangenome_node;
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
        registerShutdownHook(graphDb);
        startTime = System.currentTimeMillis();
        genomeDb = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH);
        indexDb = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH, "sorted");
        K_SIZE = indexDb.get_K();
        genome_scanner = new SequenceScanner(genomeDb, 1, 1, K_SIZE, indexDb.get_pre_len());
        address = new int[4];
        seq = new StringBuilder();
        try (Transaction tx = graphDb.beginTx()) {
            pangenome_node = graphDb.findNodes(pangenome_label).next();
            if (pangenome_node == null) {
                System.out.println("Can not locate database node!");
                System.exit(1);
            }
            K_SIZE = (int) pangenome_node.getProperty("k_mer_size");
            try {
                in = new BufferedReader(new FileReader(PATH_TO_THE_GENOME_NUMBERS_FILE));
                while (in.ready()) {
                    genome_number = in.readLine().trim();
                    if (genome_number.equals(""))
                        continue;
                    try{
                        address[0] = Integer.parseInt(genome_number);
                    }catch(NumberFormatException e){
                        System.out.println(genome_number + "is not a valid genome number.");
                        continue;
                    }
                    if (address[0] < 1 || address[0] > genomeDb.num_genomes){
                        System.out.println(genome_number + "is not a valid genome number.");
                        continue;
                    }
                    System.out.println("Reconstructing genome " + genome_number + "...");
                    try {
                        out = new BufferedWriter(new FileWriter(PATH_TO_THE_PANGENOME_DATABASE + "/Genome_" + genome_number + ".fasta"));
                        for (address[1] = 1; address[1] <= genomeDb.num_sequences[address[0]]; ++address[1]) {
                            System.out.println("Sequence " + address[1] + " length = " + genomeDb.sequence_length[address[0]][address[1]]);
                            //address[2] = 1;
                            //address[3] = (int)genomeDb.sequence_length[address[0]][address[1]];
                            //start = locate(address, K);
                            //extract_sequence(seq, start, address, K);
                            out.write(">" + genomeDb.sequence_titles[address[0]][address[1]] + "\n");
                            address[2] = 0;
                            address[3] = (int)genomeDb.sequence_length[address[0]][address[1]] - 1;
                            seq.setLength(0);
                            genome_scanner.get_sub_sequence(seq, address, true);
                            write_fasta(out, seq.toString(), 80);
                            seq.setLength(0);
                        }
                        out.close();
                    } catch (IOException e) {
                        System.out.println(e.getMessage());
                        System.exit(1);
                    }
                }
                in.close();
            } catch (IOException ioe) {
                System.out.println("Failed to read file names!");
                System.exit(1);
            }
            tx.success();
        }
        System.out.println("Genomes were stored in the database directory.");
        graphDb.shutdown();
        genomeDb.close();
    }

    public void retrieve_synteny(String genome) {
        if (! new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH).exists()) {
            System.out.println("No database found in " + PATH_TO_THE_PANGENOME_DATABASE);
            System.exit(1);
        }
        int g;
        Node pangenome_node;
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
        registerShutdownHook(graphDb);
        startTime = System.currentTimeMillis();
        genomeDb = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH);
        try (Transaction tx = graphDb.beginTx()) {
            pangenome_node = graphDb.findNodes(pangenome_label).next();
            if (pangenome_node == null) {
                System.out.println("Can not locate database node!");
                System.exit(1);
            }
            K_SIZE = (int) pangenome_node.getProperty("k_mer_size");
            tx.success();
        }
        //try {
            g = Integer.parseInt(genome.trim());
            System.out.println("Reconstructing synteny map between genome " + genome + " and the rest :");
            extract_synteny(g);
        /*}catch(NumberFormatException e){
            System.out.println("Invalid genome number!");
            System.exit(1);
        }*/
        System.out.println("Synteny files are ready in " + PATH_TO_THE_PANGENOME_DATABASE);
        graphDb.shutdown();
        genomeDb.close();
    }
    
    public void extract_synteny(int genome1) {
        Relationship rel;
        Node node, neighbor;
        IndexPointer start;
        try{
            BufferedWriter[] out_fwd = new BufferedWriter[genomeDb.num_genomes + 1];
            BufferedWriter[] out_rev = new BufferedWriter[genomeDb.num_genomes + 1];
            String formatStr = "%8s%10s%10s\n";
            int[] loc2;
            int genome2, loc1, seq1, seq2, i, node_len;
            long offset1, offset2;
            char side1, side2;
            String origin1;
            for (seq1 = 1; seq1 <= genomeDb.num_sequences[genome1]; ++seq1) {
                try (Transaction tx = graphDb.beginTx()) {
                    loc1 = 0;
                    start = locate(graphDb, new int[]{genome1, seq1, loc1 + 1}, K_SIZE);
                    origin1 = "G" + genome1 + "S" + seq1;
                    node = graphDb.getNodeById(start.node_id);
                    side1 = start.canonical ? 'F' : 'R';
                    while (true) {
                        //System.out.println(loc1);
                        node_len = (int) node.getProperty("length");
                        for (Relationship r: node.getRelationships(Direction.INCOMING)){
                            side2 = r.getType().name().charAt(1);
                            for (String origin2: r.getPropertyKeys()){
                                loc2 = (int[])r.getProperty(origin2);
                                genome2 = Integer.parseInt(origin2.split("S")[0].substring(1));
                                seq2 = Integer.parseInt(origin2.split("S")[1]);
                                offset2 = genomeDb.sequence_offset[genome2][seq2];
                                offset1 = genomeDb.sequence_offset[genome1][seq1];
                                if (side1 == side2){
                                    if (out_fwd[genome2] == null){
                                        out_fwd[genome2] = new BufferedWriter(new FileWriter(PATH_TO_THE_PANGENOME_DATABASE + "/F_"+ genome1 + "_" + genome2 + ".smf"));
                                        out_fwd[genome2].write("> " + genomeDb.genome_names[genome2].split("\\s")[0] + "\n");
                                    }
                                    for (i = 0; i < loc2.length; ++i){
                                        out_fwd[genome2].write(String.format(formatStr, offset1 + loc1 + 1, offset2 + loc2[i] + 1, node_len));
                                        }
                                } else {
                                    if (out_rev[genome2] == null){
                                        out_rev[genome2] = new BufferedWriter(new FileWriter(PATH_TO_THE_PANGENOME_DATABASE + "/R_"+ genome1 + "_" + genome2 + ".smf"));
                                        out_rev[genome2].write("> " + genomeDb.genome_names[genome2].split("\\s")[0] + " Reverse\n");
                                    }
                                    for (i = 0; i < loc2.length; ++i)
                                        out_rev[genome2].write(String.format(formatStr, offset1 + loc1 + (side1 == 'F' ? 1 : node_len), offset2 + loc2[i] + (side2 == 'F' ? 1 : node_len), node_len));
                                }
                            }
                        }
                        loc1 += node_len - K_SIZE + 1;
                        rel = get_outgoing_edge(node, origin1, loc1);
                        if (rel == null)
                            break;
                        else
                            neighbor = rel.getEndNode();
                        node = neighbor;
                        side1 = rel.getType().name().charAt(1);
                    } // while
                    System.out.print("\rSequence " + seq1 + " / " + genomeDb.num_sequences[genome1] + " finished.");
                    tx.success();
                }
            }
            System.out.println();
            for (genome2 = 1; genome2 <= genomeDb.num_genomes; ++genome2){
                if (out_fwd[genome2] != null)
                    out_fwd[genome2].close();
                if (out_rev[genome2] != null)
                    out_rev[genome2].close();  
            }
        } catch (IOException ioe) {
            System.out.println("Failed to read file names!");
            System.exit(1);
        }
    }
    
    /**
     * Extracts the genomic region belonging to the specified sequence starting at th specified node.
     * 
     * @param seq Will contains the sequence after function ends.
     * @param start_ptr A pangenome pointer which points to the node where the sequence starts.
     * @param address An array determining {genome, sequence, begin, end}properties of the sequence.
     */
    public void extract_sequence(StringBuilder seq, IndexPointer start_ptr, int[] address) {
        Relationship rel;
        Node neighbor, node;
        int[] addr = new int[]{address[0],address[1],address[2],address[3]};
        int begin = addr[2] - 1, end = addr[3] - 1;
        int len = 0, node_len, neighbor_len, seq_len, position;
        String rel_name, origin = "G" + address[0] + "S" + address[1];
        seq_len = end - begin + 1;
        seq.setLength(0);
        position = start_ptr.offset;
        node = graphDb.getNodeById(start_ptr.node_id);
        node_len = (int) node.getProperty("length");
    // Takes the part of the region lies in the first node of the path that region takes in the graph    
        if (start_ptr.canonical) {
            if (position + seq_len - 1 <= node_len - 1) { // The whole sequence lies in this node
                len += append_fwd(seq, (String) node.getProperty("sequence"), position, position + seq_len - 1);
            } else {
                len += append_fwd(seq, (String) node.getProperty("sequence"), position, node_len - 1);
            }
        } else {
            if (position - (seq_len - 1) >= 0) { // The whole sequence lies in this node
                len += append_rev(seq, (String) node.getProperty("sequence"), position - (seq_len - 1), position);
            } else {
                len += append_rev(seq, (String) node.getProperty("sequence"), 0, position);
            }
        }
    //  traverse the path of the region   
        while (len < seq_len) {
            //System.out.println(node.getId()+" "+len + " " + seq_len);
            addr[2] = (begin + len) - K_SIZE + 1;
            rel = get_outgoing_edge(node, origin, addr[2]);
            neighbor = rel.getEndNode();
            rel_name = rel.getType().name();
            neighbor_len = (int) neighbor.getProperty("length");
            if (rel_name.charAt(1) == 'F') {// Enterring forward side
                if (len + neighbor_len - K_SIZE + 1 > seq_len) // neighbor is the last node of the path
                    len += append_fwd(seq, (String) neighbor.getProperty("sequence"), K_SIZE - 1, seq_len - len + K_SIZE - 2);
                else 
                    len += append_fwd(seq, (String) neighbor.getProperty("sequence"), K_SIZE - 1, neighbor_len - 1);
            }else{ // Enterring reverse side
                if (len + neighbor_len - K_SIZE + 1 > seq_len) // neighbor is the last node of the pat
                    len += append_rev(seq, (String) neighbor.getProperty("sequence"), neighbor_len - K_SIZE - (seq_len - len) + 1, neighbor_len - K_SIZE);
                else 
                    len += append_rev(seq, (String) neighbor.getProperty("sequence"), 0, neighbor_len - K_SIZE);
            }
            node = neighbor;
        } // while
    }
  
    /**
     * Give the next node of the path to be traversed through.
     * 
     * @param current_node The current node of the path we are located at. 
     * @param origin 
     * @param address An array which determine the genome, sequence and position of the desirable outgoing edge.
     * @return The outgoing edge.
     */
    public static Relationship get_outgoing_edge(Node current_node, String origin, int pos) {
        int[] occurrence;
        for (Relationship r_out : current_node.getRelationships(Direction.OUTGOING)) {
            occurrence = (int[])r_out.getProperty(origin, null);
            if (occurrence != null) {
                if (Arrays.binarySearch(occurrence, pos) >= 0)
                    return r_out;
            }
        }
        return null;
    }
    
    /**
     * Appends substring s[from..to] to the string builder.
     * @param seq String builder.
     * @param s The string.
     * @param from Start position of the substring.
     * @param to Stop position of the substring.
     * @return The length of substring appended to the string builder.
     */
    public static int append_fwd(StringBuilder seq, String s, int from, int to) {
        for (int i = from; i <= to; ++i)
            seq.append(s.charAt(i));
        return to - from + 1;
    }
    
    /**
     * Appends the reverse complement of substring s[from..to] to the string builder.
     * @param seq String builder.
     * @param s The string.
     * @param from Start position of the substring.
     * @param to Stop position of the substring.
     * @return The length of substring appended to the string builder.
     */
    public static int append_rev(StringBuilder seq, String s, int from, int to) {
        for (int i = to; i >= from; --i)
            seq.append(complement(s.charAt(i)));
        return to - from + 1;
    }
    
    /**
     * Returns a pangenome pointer pointing to the specified genomic position.
     * 
     * @param address An integer array lile {genome_number, sequence_number, begin_position, end_position}
     * @return A pointer to the genomic position in the pangenome
     */
    public static IndexPointer locate(GraphDatabaseService graphDb, int[] addr, int K_SIZE) {
        int node_start_pos, low, high, mid , node_len, genomic_pos;
        boolean forward;
        Node node, neighbor, seq_node;
        Relationship rel;
        String anchor_sides, origin = "G" + addr[0] + "S" + addr[1];
        long[] anchor_nodes;
        int[] anchor_positions;
        int[] address = Arrays.copyOf(addr,addr.length);
        genomic_pos = address[2] - 1;
        seq_node = graphDb.findNode(sequence_label, "identifier", address[0]+"_"+address[1]);
        anchor_nodes = (long[]) seq_node.getProperty("anchor_nodes");
        anchor_positions = (int[]) seq_node.getProperty("anchor_positions");
        anchor_sides = (String) seq_node.getProperty("anchor_sides");
    // Find the immediate preceding anchor_node, searching in the sorted array of anchor positions.      
        for (low = 0, high = anchor_sides.length() - 1, mid = (low + high) / 2; low <= high; mid = (low + high) / 2) {
            if (genomic_pos < anchor_positions[mid]) {
                high = mid - 1;
            } else if (genomic_pos > anchor_positions[mid]) {
                low = mid + 1;
            } else {
                break;
            }
        }
        if (genomic_pos < anchor_positions[mid]) {
            --mid;
        }
        forward = anchor_sides.charAt(mid) == 'F';
        try (Transaction tx = graphDb.beginTx()) {
            node = graphDb.getNodeById(anchor_nodes[mid]);
            node_start_pos = anchor_positions[mid];
            node_len = (int) node.getProperty("length");
        // Traverse the pangenome from the anchor node until reach to the target
            while (node_start_pos + node_len <= genomic_pos) 
            {
                address[2] = node_start_pos + node_len - K_SIZE + 1;
                rel = get_outgoing_edge(node, origin, address[2]);
                if (rel == null){
                    System.out.println("Failed to locate address : " + address[0] + " " + address[1] + " "+ address[2]);
                    break;
                }
                neighbor = rel.getEndNode();
                forward = rel.getType().name().charAt(1) == 'F';
                node_start_pos += node_len - K_SIZE + 1;
                node = neighbor;
                node_len = (int) node.getProperty("length");
            }
            tx.success();
        }
        return new IndexPointer(node.getId(), forward, forward ? genomic_pos - node_start_pos : node_len - 1 - (genomic_pos - node_start_pos), -1l);
    }
    
    /***
     * Creates an edge between source and destination nodes.
     * 
     * @param src Source node
     * @param des Destination node
     * @param edge_type One of the four possible edge types: FF, FR, RF, RR
     * @param address Specifies which genomic address the edge points to. 
     * @return The newly created edge
     */
    private void connect(Node src, Node des, RelationshipType edge_type) {
        if (DEBUG) System.out.println("connect "+src.getId()+" "+edge_type.name()+" "+des.getId());
        src.createRelationshipTo(des, edge_type);
    }
    
    /**
     * Splits a node at a specified position by creating a new node called split_node as a part separated from the node.
     * @param node The node which should be split.
     * @param pos The position of the split with respect to the start on the node.
     * @return The newly created split node.
     */
    private Node split(Node node, int pos) {
        int split_len, node_len;
        int i, s_id, gen, seq, loc,starts_at;
        long inx,split_first_kmer,node_last_kmer=0;
        int[] address;
        Node neighbor, split_node;
        Relationship rel;
        address = (int[]) node.getProperty("address");
        gen = address[0];
        seq = address[1];
        loc = address[2];
        node_len = (int) node.getProperty("length");
        ++num_nodes;
        split_node = graphDb.createNode(nucleotide_label);
        address[0] = gen;
        address[1] = seq;
        address[2] = loc + pos;
        split_node.setProperty("address", address);
        split_len = node_len - pos;
        split_node.setProperty("length", split_len);
    // Updates the edges comming from gene level to the node.    
        for (Relationship r : node.getRelationships(RelTypes.starts, Direction.INCOMING)) {
            starts_at = (int)r.getProperty("offset");
            if (starts_at >= pos) {
                rel = r.getStartNode().createRelationshipTo(split_node, RelTypes.starts);
                rel.setProperty("offset", starts_at - pos);
                rel.setProperty("forward", r.getProperty("forward"));
                rel.setProperty("genomic_position", r.getProperty("genomic_position"));
            r.delete();
            } 
        }        
        // Updating the Kmers chain in the index  
        node_last_kmer=indexDb.find(genome_scanner.make_kmer(gen,seq,loc+pos-1));
        split_first_kmer=indexDb.get_next_index(node_last_kmer);
        indexDb.put_next_index(-1L, node_last_kmer); 
        split_node.setProperty("first_kmer",split_first_kmer);
        split_node.setProperty("last_kmer",node.getProperty("last_kmer"));
        s_id=(int)split_node.getId();
        for(i=0,inx=split_first_kmer;inx!=-1L;inx=indexDb.get_next_index(inx),++i) // update kmer coordinates
        {
            indexDb.put_node_id(s_id, inx);
            indexDb.put_position(i, inx);
        }  
        // Moving forward-outgoing and reverse-incoming edges from node to split node.    
        for (Relationship r : node.getRelationships(Direction.OUTGOING,RelTypes.FR,RelTypes.FF)) {
            neighbor = r.getEndNode();
            if (neighbor.equals(node)) 
                neighbor = r.isType(RelTypes.FF) ? node : split_node;
            connect(split_node, neighbor, r.getType());
            r.delete();
        }
        for (Relationship r : node.getRelationships(Direction.INCOMING,RelTypes.RR,RelTypes.FR)) {
            neighbor = r.getStartNode();
            if (neighbor.equals(node)) 
                neighbor = r.isType(RelTypes.RR) ? node : split_node;
            connect(neighbor, split_node, r.getType());
            r.delete();
        }
    //  Connecting node to split node
        if (node.hasRelationship(Direction.INCOMING, RelTypes.FF, RelTypes.RF)){
            connect(node, split_node, RelTypes.FF);
            ++num_edges;
        }
        if (split_node.hasRelationship(Direction.INCOMING, RelTypes.FR, RelTypes.RR)){
            connect(split_node ,node, RelTypes.RR);
            ++num_edges;
        }
        node.setProperty("last_kmer",node_last_kmer);
        node.setProperty("length", pos + K_SIZE - 1);
        return split_node;
    }
    
    /**
     * Creates and extends a new node till reach to a previously visited K-mer or a degenerate region.
     * 
     * @param node The node to be extended.
     */
    private void create_extend() {
        int[] address;
        long node_id,last_kmer;
        address = genome_scanner.get_address();
        address[2] -= K_SIZE - 1;
        ++num_nodes;
        Node new_node = graphDb.createNode(nucleotide_label);
        node_id = new_node.getId();
        if (DEBUG) System.out.println("create "+new_node.getId());
        new_node.setProperty("address", address);
        new_node.setProperty("length", K_SIZE);
        new_node.setProperty("last_kmer",genome_scanner.get_curr_index());
        new_node.setProperty("first_kmer",genome_scanner.get_curr_index());
    // Set the pointer to the Kmer in the pointer database    
        indexDb.put_pointer(new_node.getId(), 0, genome_scanner.get_curr_kmer().get_canonical(), -1l, genome_scanner.get_curr_index());
        connect(curr_node ,new_node, RelTypes.values()[curr_side*2]);
        ++num_edges;
        curr_node = new_node;
        curr_side = 0;

        if (DEBUG) System.out.println("extending node "+curr_node.getId());
        int begin, len;
        last_kmer=(long)curr_node.getProperty("last_kmer");
        boolean broke;
        Node degenerate_node = null;
        len = (int) curr_node.getProperty("length");
        broke = false;
        while (genome_scanner.get_position() < genome_scanner.get_sequence_length() - 1) { // Not reached to the end of the sequence
            if (DEBUG) System.out.println("extend " + genome_scanner.get_position());
            if (genome_scanner.get_code(1) > 3) { // hit a degenerate region
                genome_scanner.next_position();
                begin = genome_scanner.get_position() - K_SIZE + 1;
                curr_node.setProperty("length", len);
                curr_node.setProperty("last_kmer",last_kmer);                
                genome_scanner.jump_forward();
                if (genome_scanner.get_position() >= genome_scanner.get_sequence_length() - 1){
                    genome_scanner.next_position();// to acheive the right length for the degenerate node    
                    finish = true;
                }                
                int[] add = genome_scanner.get_address();
                add[2] = begin;
                degenerate_node = create_degenerate(add);
                connect(curr_node ,degenerate_node, RelTypes.FF); 
                ++num_edges;
                curr_node = degenerate_node;
                break;
            } else {
                genome_scanner.next_position();
                genome_scanner.get_curr_kmer().next_kmer(genome_scanner.get_code(0));
                if (SHOW_KMERS) System.out.println(genome_scanner.get_curr_kmer().toString());
            }
            genome_scanner.set_curr_index(indexDb.find(genome_scanner.get_curr_kmer()));
            if (indexDb.get_node_id(genome_scanner.get_curr_index()) == -1L) {
                indexDb.put_next_index(genome_scanner.get_curr_index(),last_kmer);
                ++len;
                indexDb.put_pointer(node_id, len - K_SIZE, genome_scanner.get_curr_kmer().get_canonical(), -1l, genome_scanner.get_curr_index());
                last_kmer=genome_scanner.get_curr_index();
            } else {
                broke = true;
                break;
            }
        }
        if (degenerate_node == null){
            new_node.setProperty("length", len);
            new_node.setProperty("last_kmer",last_kmer);
        }
        if (!broke && genome_scanner.get_position() >= genome_scanner.get_sequence_length() - 1) {// because the last kmer is somewhere in the graph we should get connected to
            finish = true;
        }
    }

    /**
     * Enters the node in which the current Kmer is found in and performs zero, one or two splits.   
     */
    private void follow_forward(IndexPointer pointer) {
        int l, pos, begin, g, s, loc, side;
        Node node, split_node1, split_node2,des, src, degenerate_node = null;
        RelationshipType rel_type;
        int[] address;
        boolean loop, repeated_edge;
        pos = pointer.offset;
        node = graphDb.getNodeById(pointer.node_id);
        if (DEBUG) System.out.println("follow_forward "+pointer.node_id+" at "+pos);
    // The first split might be done to seperate the part we need to enter in.
        if (pos > 0) {  
            if (DEBUG) System.out.println("first_split "+node.getId()+" at "+pos);
            split_node1 = split(node, pos);
            if (loop = (curr_node.equals(node) && curr_side == 0)) 
                src = split_node1;
            else 
                src = curr_node;
            node = split_node1; // Note : assigning reference variables is dangerous! if split_node changes node will change as well.
        } else {
            split_node1 = node;
            if (loop = (curr_node.equals(node) && curr_side == 0)) 
                src = split_node1;
            else 
                src = curr_node;
        }
        des = split_node1;
        side=curr_side*2;
        curr_side = 0;
        l = (int) node.getProperty("length") - K_SIZE;
        address = (int[]) node.getProperty("address");
        g = address[0];
        s = address[1];
        loc = address[2];
    // Follow the shared part
        for (pos = 0; pos <= l && genome_scanner.get_position() <= genome_scanner.get_sequence_length() - 1 && genome_scanner.get_code(g, s, loc + pos + K_SIZE - 1) == genome_scanner.get_code(0); ++pos) {
            genome_scanner.next_position();
         // If hit a degenarate region aplit and branch to a degenerate node 
            if (genome_scanner.get_position() <= genome_scanner.get_sequence_length() - 1 && genome_scanner.get_code(0) > 3) {
                begin = genome_scanner.get_position() - K_SIZE + 1;
                genome_scanner.jump_forward();
                if (genome_scanner.get_position() >= genome_scanner.get_sequence_length() - 1){
                    genome_scanner.next_position();// to acheive the right length for the degenerate node    
                    finish = true;
                }                
                if (pos + 1 <= l) {
                    split_node2 = split(node, pos + 1);
                    if (loop)
                        src = split_node2;
                }
                int[] add = genome_scanner.get_address();
                add[2] = begin;
                degenerate_node = create_degenerate(add);
                connect(node ,degenerate_node, RelTypes.FF);
                ++num_edges;
                break;
            }
        }
        if (genome_scanner.get_position() == genome_scanner.get_sequence_length()) {
            finish = true;
        } else if (degenerate_node == null) // build the Kmer of difference 
            initialize(genome_scanner.get_position() - K_SIZE + 1);
    //  A second split might be needed   
        if (degenerate_node == null && pos <= l) {
            if (DEBUG) System.out.println("second_split "+node.getId()+" at "+pos);
            split_node2 = split(node, pos);
            if (loop)
                src = split_node2;
        }
    // connect the current node before doing splits to the split_node1    
        rel_type = RelTypes.values()[side];
        repeated_edge = false;
        for (Relationship r: src.getRelationships(rel_type, Direction.OUTGOING))
            if (r.getEndNode().equals(des)){
                repeated_edge = true;
                break;
            }
        if (!repeated_edge){
            connect(src ,des, rel_type);
            ++num_edges;
        }
        if (degenerate_node != null) {
            curr_side = 0; // not really needed
            curr_node = degenerate_node;
        } else
            curr_node = node;
    }
    
    /**
     * Enters the forward side of node in which the current Kmer found and performs zero, one or two splits.   
     */
    private void follow_reverse(IndexPointer pointer) {
        int pos, begin, g, s, loc, side;
        int[] address;
        Node node, split_node1, split_node2 ,des, src, degenerate_node = null;
        boolean loop, first_split = false, repeated_edge;
        pos = pointer.offset;
        node = graphDb.getNodeById(pointer.node_id);
        if (DEBUG) System.out.println("follow_reverse "+pointer.node_id+" at "+pos);
        RelationshipType rel_type;
        split_node2 = node; //if the second split does not happens remains unchanged
        if (pos < (int) node.getProperty("length") - K_SIZE) {
            if (DEBUG) System.out.println("first_split "+node.getId()+" at "+(pos+1));
            first_split = true;
            split_node1 = split(node, pos+1);
            if (loop = curr_node.equals(node) && curr_side == 0) // might be in reverse side due to a follow reverse
                src = split_node1;
            else 
                src = curr_node;
        } else {
            split_node1 = node;
            if (loop = curr_node.equals(node) && curr_side == 0)
                src = split_node1;
            else 
                src = curr_node;
        }
        des = node;
        side=curr_side*2+1;
        curr_side = 1;
        address = (int[]) node.getProperty("address");
        g = address[0];
        s = address[1];
        loc = address[2];
        for (pos = (int) node.getProperty("length") - K_SIZE; pos >= 0 && genome_scanner.get_position() <= genome_scanner.get_sequence_length() - 1 && genome_scanner.get_code(g, s, loc + pos) == genome_scanner.get_complement_current_code(0); --pos) {
            genome_scanner.next_position();
            if (genome_scanner.get_position() <= genome_scanner.get_sequence_length() - 1 && genome_scanner.get_code(0) > 3) {
                begin = genome_scanner.get_position() - K_SIZE + 1;
                genome_scanner.jump_forward();
                if (genome_scanner.get_position() >= genome_scanner.get_sequence_length() - 1){
                    genome_scanner.next_position();// to acheive the right length for the degenerate node    
                    finish = true;
                }                
                if (pos > 0) {
                    split_node2 = split(node, pos);
                    des = split_node2;
                    if (!first_split && loop)
                        src = split_node2;
                }
                int[] add = genome_scanner.get_address();
                add[2] = begin;
                degenerate_node = create_degenerate(add);
                connect(split_node2, degenerate_node, RelTypes.RF);
                ++num_edges;
                break;
            }
        }
        if (genome_scanner.get_position() == genome_scanner.get_sequence_length()) {
            finish = true;
        } else if (degenerate_node == null) // build the Kmer of difference 
            initialize(genome_scanner.get_position() - K_SIZE + 1);
        if (degenerate_node == null && pos >= 0) {
            if (DEBUG) System.out.println("second_split "+node.getId()+" at "+ (pos + 1));
            split_node2 = split(node, pos+1);
            des = split_node2;
            if (!first_split && loop)
                src = split_node2;
        }
        rel_type = RelTypes.values()[side];
        repeated_edge = false;
        for (Relationship r: src.getRelationships(rel_type, Direction.OUTGOING))
            if (r.getEndNode().equals(des)){
                repeated_edge = true;
                break;
            }
        if (!repeated_edge){
            connect(src ,des, rel_type);
            ++num_edges;
        }
        if (degenerate_node != null) {
            curr_side = 0;
            curr_node = degenerate_node;
        } else
            curr_node = split_node2;
    }

    /**
     * creates a degenerate node starting at "begin" ending at position-1.
     * @param address The genomic position of the region
     */
    private Node create_degenerate(int[] address) {
        ++num_degenerates;
        ++num_nodes;
        Node degenerate_node = graphDb.createNode(degenerate_label);
        degenerate_node.addLabel(nucleotide_label);
        if (DEBUG) System.out.println("create_degenerate:"+degenerate_node.getId()+" position:"+genome_scanner.get_position()+" begin:"+address[2]);
        degenerate_node.setProperty("address", address);
        degenerate_node.setProperty("length", genome_scanner.get_position() - address[2]);
        return degenerate_node;
    }
    
    private void initialize(int start){
        Node degenerate_node;
        if (genome_scanner.initialize_kmer(start) < 0){// start with a degenerate kmer
            if (genome_scanner.get_position() >= genome_scanner.get_sequence_length() - 1){
                genome_scanner.next_position();// to acheive the right length for the degenerate node    
                finish = true;
            }                
            int[] add = genome_scanner.get_address();
            add[2] = 0;
            degenerate_node = create_degenerate(add);
            connect(curr_node ,degenerate_node, RelTypes.values()[curr_side*2]);
            ++num_edges;
            curr_node = degenerate_node;  
        }
    }
    
    /**
     * 
     Constructs the pangenome out of the provided input sequences.
     */
    void construct_pangenome(Node pangenome_node) {
        int trsc = 0;
        Node genome_node, sequence_node;
        IndexPointer pointer = new IndexPointer();
        phaseTime = System.currentTimeMillis();
        Transaction tx = graphDb.beginTx();
        try {
            while (!genome_scanner.end_of_scan()) {
                System.out.println("Processing genome " + genome_scanner.get_genome() + " (" + genomeDb.genome_names[genome_scanner.get_genome()] +") :");
                genome_node = graphDb.createNode(genome_label);
                genome_node.setProperty("path", genomeDb.genome_names[genome_scanner.get_genome()]);
                genome_node.setProperty("number", genome_scanner.get_genome());
                genome_node.setProperty("num_sequences", genomeDb.num_sequences[genome_scanner.get_genome()]);
                genome_node.setProperty("date", new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(new Date()));
                pangenome_node.createRelationshipTo(genome_node, RelTypes.has);
                while (!genome_scanner.end_of_genome()) {
                    sequence_node = graphDb.createNode(sequence_label);
                    sequence_node.setProperty("genome", genome_scanner.get_genome());
                    sequence_node.setProperty("number", genome_scanner.get_sequence());
                    sequence_node.setProperty("identifier", genome_scanner.get_genome() + "_" + genome_scanner.get_sequence());
                    sequence_node.setProperty("title", genome_scanner.get_title());
                    sequence_node.setProperty("length", genome_scanner.get_sequence_length());
                    sequence_node.setProperty("offset", genome_scanner.get_offset());
                    genome_node.createRelationshipTo(sequence_node, RelTypes.has);
                    finish = false;
                    System.out.println("sequence " + genome_scanner.get_sequence() + "/" + genomeDb.num_sequences[genome_scanner.get_genome()] + 
                    " of genome " + genome_scanner.get_genome() + "\tlength=" + genome_scanner.get_sequence_length());
                    curr_node = sequence_node;
                    curr_side = 0;
                    if (genome_scanner.get_sequence_length() >= K_SIZE){
                        initialize(0);
                        while (!finish) {
                            genome_scanner.set_curr_index(indexDb.find(genome_scanner.get_curr_kmer()));
                            indexDb.get_pointer(pointer, genome_scanner.get_curr_index());
                            if (pointer.node_id == -1L) // kmer is new
                                create_extend();
                            else if (genome_scanner.get_curr_kmer().get_canonical() ^ pointer.canonical)// if sides don't agree
                                follow_reverse(pointer);
                            else 
                                follow_forward(pointer);
                            ++trsc;
                            if (trsc >= MAX_TRANSACTION_SIZE){    
                                tx.success();
                                tx.close();
                                tx = graphDb.beginTx();
                                trsc = 0;
                            }
                        }
                        connect(curr_node, sequence_node, RelTypes.values()[curr_side*2]);// to point to the last k-mer of the sequence located in the other strand
                        ++num_edges;
                    }
                    genome_scanner.next_sequence();
                }//sequences
                System.out.println((System.currentTimeMillis() - phaseTime) / 1000 + " seconds elapsed.");
                genome_scanner.next_genome();
            }//genomes
            tx.success();
        } finally {
            tx.close();
        }
        add_sequence_properties();
        localize_nodes();
    }
    
    /**
     * To add list of anchor nodes, anchor sides and anchor positions to each sequence_node.
     * These are used for locating genomic regions very quickly. 
     */
    void localize_nodes() {
        ResourceIterator<Node> sequence_iterator;
        LinkedList<Node> sequence_nodes;
        int trsc = 0, i, len, m, neighbor_length = 0;
        char node_side, neighbor_side;
        long length, distance;
        long[] anchor_nodes;
        int[] anchor_positions;
        int[] initial_coordinate = new int[1];
        Node node, neighbor = null, sequence_node;
        StringBuilder nds = new StringBuilder();
        StringBuilder pos = new StringBuilder();
        StringBuilder sds = new StringBuilder();
        String[] ids_list, posis_list;
        String rel_name,origin;
        int[] positions;
        int[] new_positions;
        int[] address = new int[3], addr = null;
        boolean is_node = false, is_degenerate = false, found = true;
        try (Transaction tx = graphDb.beginTx()){
            sequence_iterator = graphDb.findNodes(sequence_label);
            sequence_nodes = new LinkedList();
            while (sequence_iterator.hasNext())
                sequence_nodes.add(sequence_iterator.next());
            sequence_iterator.close();
            tx.success();
        }

        Transaction tx = graphDb.beginTx();
        try {
            while (!sequence_nodes.isEmpty()){
                sequence_node = sequence_nodes.remove();
                origin = "G" + ((String)sequence_node.getProperty("identifier")).replace('_', 'S');
                address[0] = (int)sequence_node.getProperty("genome");
                address[1] = (int)sequence_node.getProperty("number");
                length = (long) sequence_node.getProperty("length");
                System.out.println("\rLocalizing sequence "+address[1] + "/" + genomeDb.num_sequences[address[0]] + " of genome " + address[0] + "                        ");
                if (length >= K_SIZE){
                    node = sequence_node;
                    node_side = 'F';
                    distance = 0;
                    for (address[2] = 0; address[2] + K_SIZE - 1 < length && found;){ // K-1 bases of the last node not added
                        //System.out.println((address[2] + K - 1)+" ? " + length);
                        found = false;
                        for (Relationship r : node.getRelationships(Direction.OUTGOING)) {
                            rel_name = r.getType().name();
                            if (rel_name.charAt(0) != node_side)
                                continue;
                            neighbor = r.getEndNode();
                            neighbor_side = rel_name.charAt(1);
                            is_node = neighbor.hasLabel(nucleotide_label) && !neighbor.hasLabel(degenerate_label);
                            is_degenerate = neighbor.hasLabel(degenerate_label);
                            if (is_node || is_degenerate){
                                addr = (int[]) neighbor.getProperty("address");
                                neighbor_length = (int) neighbor.getProperty("length");
                            }
                            //System.out.println(node.getId()+" "+address[2]+" "+node_side);
                            //System.out.println(neighbor.getId()+" "+addr[2]+" "+neighbor_side);
                            if ((is_node && genome_scanner.compare(address, addr, K_SIZE - 1,
                                    neighbor_side == 'F' ? K_SIZE - 1 : neighbor_length - K_SIZE, 1, neighbor_side == 'F'))
                                    || (is_degenerate && Arrays.equals(addr, address))) {
                                //System.out.println("found "+address[2]+" "+neighbor.getId());
                                found = true;
                                positions = (int[]) r.getProperty(origin, null);
                                if (positions != null) {
                                    len = positions.length;
                                    new_positions = new int[len + 1];
                                    for (i = 0; i < len; ++i) {
                                        new_positions[i] = positions[i];
                                    }
                                    new_positions[i] = address[2];
                                    r.setProperty(origin, new_positions);
                                    if (len > 1000)
                                        neighbor.addLabel(low_complexity_label);
                                } else {
                                    initial_coordinate[0] = address[2];
                                    r.setProperty(origin, initial_coordinate);
                                }
                                if (address[2] >= distance) {
                                    nds.append(neighbor.getId()).append(" ");
                                    sds.append(neighbor_side);
                                    pos.append(address[2]).append(" ");
                                    distance += ANCHORS_DISTANCE;
                                }
                                address[2] = address[2] + neighbor_length - K_SIZE + 1;
                                node = neighbor;
                                node_side = neighbor_side;
                                break;
                            }
                        }
                        ++trsc;
                        if (trsc >= 1000 * MAX_TRANSACTION_SIZE){
                            tx.success();
                            tx.close();
                            tx = graphDb.beginTx();
                            trsc = 0;
                        }
                        //System.out.print("%" + address[2] * 100 / length + "\t\r");
                    }
                    if (!found) {
                        System.out.println("Could not locate position " + address[2] + " from node ID=" + node.getId());
                        System.exit(1);
                    }
                    m = sds.length();
                    ids_list = nds.toString().split("\\s");
                    posis_list = pos.toString().split("\\s");
                    anchor_nodes = new long[m];
                    anchor_positions = new int[m];
                    for (i = 0; i < m; ++i) {
                        anchor_nodes[i] = Long.valueOf(ids_list[i]);
                        anchor_positions[i] = Integer.valueOf(posis_list[i]);
                    }
                    sequence_node.setProperty("anchor_nodes", anchor_nodes);
                    sequence_node.setProperty("anchor_positions", anchor_positions);
                    sequence_node.setProperty("anchor_sides", sds.toString());
                    nds.setLength(0);
                    pos.setLength(0);
                    sds.setLength(0);
                }
            }//while
            System.out.println((System.currentTimeMillis() - phaseTime) / 1000 + " seconds elapsed.");
            tx.success();
        } finally {
           tx.close();
        }
        System.out.println();
    }

    /**
     * Extracts the sequence of the nodes from the genome database and store it as a property of the node.
     */
    void add_sequence_properties() {
        int trsc = 0, node_length;
        int[] addr;
        Node node;
        ResourceIterator<Node> nodes_iterator;
        LinkedList<Node> nodes = new LinkedList();
        StringBuilder sequence = new StringBuilder();
        System.out.println("Adding sequence to the nodes...");
        try(Transaction tx = graphDb.beginTx()){
            nodes_iterator = graphDb.findNodes(nucleotide_label);
            while (nodes_iterator.hasNext()){
                nodes.add(nodes_iterator.next());
            }
            tx.success();
            nodes_iterator.close();
        }
        Transaction tx = graphDb.beginTx();
        try {
        //num_bases = K - 1; // for the missed overlapped of the last node of each sequence which will not be stored 
            while (!nodes.isEmpty()){
                node = nodes.remove();
                addr = (int[]) node.getProperty("address");
                node_length = (int) node.getProperty("length");
                num_bases += node_length;
                sequence.setLength(0);
                //num_bases += node_length - K + 1;
                //node.setProperty("sequence", genomeDb.get_sequence(addr[0], addr[1], addr[2], node_length - K + 1, true).toString());
                genome_scanner.get_sub_sequence(sequence, addr[0], addr[1], addr[2], node_length, true);
                node.setProperty("sequence", sequence.toString());
                ++trsc;
                if (trsc >= MAX_TRANSACTION_SIZE){
                    tx.success();
                    tx.close();
                    tx = graphDb.beginTx();
                    trsc = 0;
                }
            }
            tx.success();
        } finally {
            tx.close();
        }
    }

    /**
     * Remove the given property from all the nodes and degenerates.
     * @param property 
     */    
    void drop_nodes_property(String property) {
        int i;
        ResourceIterator<Node> nodes_iterator;
        LinkedList<Node> nodes = new LinkedList();
        Node node;
        try(Transaction tx = graphDb.beginTx()){
            nodes_iterator = graphDb.findNodes(nucleotide_label);
            while (nodes_iterator.hasNext()){
                nodes.add(nodes_iterator.next());
            }
            tx.success();
            nodes_iterator.close();
        }
        while (!nodes.isEmpty()) {
            try (Transaction tx = graphDb.beginTx()) {
                for (i = 0; i < MAX_TRANSACTION_SIZE && !nodes.isEmpty(); ++i) {
                    node = nodes.remove();
                    node.removeProperty(property);
                }
                tx.success();
            }
        }
    }

    /**
     * Remove the occurrence arrays of the edges.
     */    
    void drop_edges_colors() {
        int i;
        ResourceIterator<Relationship> rels;
        Relationship r;
        try (Transaction tx = graphDb.beginTx()) {
            rels = graphDb.getAllRelationships().iterator();
            tx.success();
        }
        while (rels.hasNext()) {
            try (Transaction tx = graphDb.beginTx()) {
                for (i = 0; i < MAX_TRANSACTION_SIZE && rels.hasNext(); ++i) {
                    r = rels.next();
                    if (r.isType(RelTypes.FF) || r.isType(RelTypes.FR) || r.isType(RelTypes.RF) || r.isType(RelTypes.RR))
                        for(String p:r.getPropertyKeys())
                            r.removeProperty(p);
                }
                tx.success();
            }
        }
        rels.close();
    }

    /**
     * Shuts down the graph database if the program halts unexpectedly.
     * 
     * @param graphDb The graph database object 
     */
    private static void registerShutdownHook(final GraphDatabaseService graphDb) {
        Runtime.getRuntime().addShutdownHook(new Thread() {
            @Override
            public void run() {
                graphDb.shutdown();
            }
        });
    }

    /**
     * Returns size of a given folder.
     * 
     * @param dir The folder File object.
     * @return Size of the folder in MegaBytes
     */
    public static long getFolderSize(File dir) {
        long size = 0;
        for (File file : dir.listFiles()) {
            if (file.isFile()) {
                // System.out.println(file.getName() + " " + file.length());
                size += file.length();
            } else {
                size += getFolderSize(file);
            }
        }
        return size / 1048576 + 1;
    }
}
