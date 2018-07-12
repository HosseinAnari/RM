/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pangenome;

import alignment.BoundedLocalSequenceAlignment;
import alignment.LocalSequenceAlignment;
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
import org.neo4j.graphdb.Lock;
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
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.neo4j.graphdb.NotFoundException;
import static pantools.Pantools.DEBUG;
import static pantools.Pantools.GAP_EXT;
import static pantools.Pantools.GAP_OPEN;
import static pantools.Pantools.INNER_READS_DISTANCE;
import static pantools.Pantools.K_SIZE;
import static pantools.Pantools.PATH_TO_THE_FIRST_SRA;
import static pantools.Pantools.PATH_TO_THE_SECOND_SRA;
import static pantools.Pantools.READS_DATABASE_PATH;
import static pantools.Pantools.OUTPUT_PATH;
import static pantools.Pantools.SHOW_KMERS;
import static pantools.Pantools.CLIP;
import static pantools.Pantools.THREADS;
import static pantools.Pantools.degenerate_label;
import static pantools.Pantools.genome_label;
import static pantools.Pantools.low_complexity_label;
import static pantools.Pantools.nucleotide_label;
import static pantools.Pantools.pangenome_label;
import static pantools.Pantools.sequence_label;

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
    private AtomicInteger[] num_genomic_mapping;
    private AtomicInteger num_mapping;
    private AtomicInteger read_number;
    //private AtomicInteger number_of_alignments;
    private GraphDatabaseService graphDb;
    private IndexDatabase indexDb;
    private SequenceDatabase genomeDb;
    private SequenceDatabase[] sequencingDb;
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
    
    public class alignment{
        public String reference;
        public hit hit;

        public alignment(String r, hit h){
            reference = r;
            hit = h;
        }
    }
    
    public class hit{
        public int genome;
        public int sequence;
        public double score;
        public int start;
        public int offset;
        public boolean forward;
        public String cigar;

        public hit(int gn, int sq, double sc, int st, int off, boolean fw, String cg){
            genome = gn;
            sequence = sq;
            score = sc;
            start = st;
            offset = off;
            forward = fw;
            cigar = cg;
        }
        
        public hit(){
        }

        public void copy(int gn, int sq, double sc, int st, int off, boolean fw, String cg){
            genome = gn;
            sequence = sq;
            score = sc;
            start = st;
            offset = off;
            forward = fw;
            cigar = cg;
        }

        public void copy(hit h){
            genome = h.genome;
            sequence = h.sequence;
            score = h.score;
            start = h.start;
            offset = h.offset;
            forward = h.forward;
            cigar = new String(h.cigar);
        }

        @Override
        public String toString(){
            return "(genome:" + genome + 
                    ",sequence:" + sequence + 
                    ",score:" + score + 
                    ",start:" + start + 
                    ",offset:" + offset + 
                    ",forward:" + forward + 
                    ",cigar:" + cigar +")";
        }
    }

    /**
     * Implements a comparator for integer arrays of size two
     */
    public class hitComparator implements Comparator<hit> {
        @Override
        public int compare(hit x, hit y) {
            if (y.score < x.score) 
                return -1;
            else if (y.score > x.score) 
                return 1;
            else if (y.cigar.length() > x.cigar.length()) 
                return -1;
            else if (y.cigar.length() < x.cigar.length()) 
                return 1;
            else if (y.sequence < x.sequence) 
                return -1;
            else if (y.sequence > x.sequence) 
                return 1;
            else if (y.start < x.start) 
                return -1;
            else if (y.start > x.start) 
                return 1;
            else
                return 0;
        }
    }    
    
    public class CountComparator implements Comparator<int[]> {
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

    public class IntComparator implements Comparator<Integer> {

        @Override
        public int compare(Integer v1, Integer v2) {
            return v1 < v2 ? -1 : v1 > v2 ? +1 : 0;
        }
    }

    public class Map implements Runnable {
        double MIN_ALIGNMENT_SCORE = 30.0;
        int KMER_SAMPLES = 20;
        int MAX_ALIGNMENT_LENGTH = 500;
        int ALIGNMENT_MODE = 0; // 0:any-best, 1:all-best
        int MINBOUND = 3;
        int MAXBOUND = 23;
        SequenceDatabase genome_database;
        IndexDatabase index_database;
        int number;
        int K;
        read[] reads;
        BufferedWriter[] mapped;
        SequenceScanner[] reads_scanner;    
        ArrayList<Integer>[][][] locations;
        IntComparator intcomp = new IntComparator();
        int[] genome_numbers;
        IndexPointer anchor_pointer;
        StringBuilder reference;
        StringBuilder sam_record;
        LinkedList<Integer>[][] sequences;
        BoundedLocalSequenceAlignment bounded_aligner;
        //LocalSequenceAlignment bounded_aligner;
        LocalSequenceAlignment aligner;
        PriorityQueue<hit>[] hits;
        hitComparator hitcomp;
        LinkedList<alignment>[] alignments;
        StringBuilder rname;
        int anchor_position;
        boolean paired_end;
        CountComparator cnt_comp;
        ArrayList<int[]> hit_counts;
        int num_exceptions = 0;
        int num_neighbors = 0;
        int num_kmers = 0;
        int num_found = 0;
        StringBuilder read_name;
        int[] address = new int[3];
        
        public Map(int n, int[] gn, BufferedWriter[] m, boolean paired) {
            int i, j, genome;
            index_database = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH);
            genome_database = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH);
            K = index_database.get_K();
            genome_scanner = new SequenceScanner(genome_database, 1, 1, K, index_database.get_pre_len());
            number = n;
            genome_numbers = Arrays.copyOf(gn, gn.length);
            mapped = m;
            paired_end = paired; 
            hitcomp = new hitComparator();
            hits = new PriorityQueue[2];
            alignments = new LinkedList[2];
            reads_scanner = new SequenceScanner[2];
            sequences = new LinkedList[2][];
            locations = new ArrayList[2][][];
            reads = new read[2];
            hits[0] = new PriorityQueue(hitcomp);
            alignments[0] = new LinkedList();
            reads_scanner[0] = new SequenceScanner(sequencingDb[0], 1, 1, K, index_database.get_pre_len());
            sequences[0] = new LinkedList[genome_database.num_genomes + 1];
            locations[0] = new ArrayList[genome_database.num_genomes + 1][];
            reads[0] = new read();
            for (i = 0; i < genome_numbers.length; ++i){
                genome = genome_numbers[i];
                locations[0][genome] = new ArrayList[genome_database.num_sequences[genome] + 1];
                for (j = 1; j <= genome_database.num_sequences[genome]; ++j)
                    locations[0][genome][j] = new ArrayList();
                sequences[0][genome] = new LinkedList();
            }
            if (paired_end){
                hits[1] = new PriorityQueue(hitcomp);
                alignments[1] = new LinkedList();
                reads_scanner[1] = new SequenceScanner(sequencingDb[1], 1, 1 ,K, index_database.get_pre_len());
                sequences[1] = new LinkedList[genome_database.num_genomes + 1];
                locations[1] = new ArrayList[genome_database.num_genomes + 1][];
                reads[1] = new read();
                for (i = 0; i < genome_numbers.length; ++i){
                    genome = genome_numbers[i];
                    locations[1][genome] = new ArrayList[genome_database.num_sequences[genome] + 1];
                    for (j = 1; j <= genome_database.num_sequences[genome]; ++j)
                        locations[1][genome][j] = new ArrayList();
                    sequences[1][genome] = new LinkedList();
                }
            }
            anchor_pointer = new IndexPointer();
            reference = new StringBuilder();
            //bounded_aligner = new LocalSequenceAlignment(GAP_OPEN, GAP_EXT, MAX_ALIGNMENT_LENGTH, CLIP, 'N');
            bounded_aligner = new BoundedLocalSequenceAlignment(GAP_OPEN, GAP_EXT, MAX_ALIGNMENT_LENGTH, MAXBOUND, CLIP, 'N');
            //bounded_aligner.initialize_bound(bound, MAX_ALIGNMENT_LENGTH);
            aligner = new LocalSequenceAlignment(GAP_OPEN, GAP_EXT, MAX_ALIGNMENT_LENGTH, CLIP, 'N');
            sam_record = new StringBuilder();
            cnt_comp = new CountComparator();
            hit_counts = new ArrayList();
            read_name = new StringBuilder();
        }

        @Override
        public void run() {
            int read_num, trsc;
            int i;
            int num_reads;
            int[] reads_to_genome;
            int genome, read_mapped = 0, counter;
            boolean mapped;
            num_reads = sequencingDb[0].num_sequences[1];
            counter = number * num_reads / 50 / THREADS;
            reads_to_genome = new int[genome_database.num_genomes + 1];
            read_num = read_number.getAndIncrement();
            while (read_num <= num_reads) {
                try(Transaction tx = graphDb.beginTx()){
                    for (trsc = 1; trsc <= 10000 && read_num <= num_reads; ++trsc, ++counter, read_num = read_number.getAndIncrement()){
                        anchor_read(read_num);
                        alignments[0].clear();
                        if (paired_end)
                            alignments[1].clear();
                        for (mapped = false, i = 0; i < genome_numbers.length; ++i){
                            genome = genome_numbers[i];
                            find_hits(genome);
                            if ((!paired_end && single_read_mapped(genome)) || (paired_end && paired_read_mapped(genome))){
                                reads_to_genome[genome]++;
                                mapped = true;
                            }
                        }
                        if (mapped)
                            ++read_mapped;
                        if (counter % (1 + num_reads / 50) == 0){
                            System.out.print("|");
                        }
                    }
                    tx.success();
                }
            }// for each read
            num_mapping.getAndAdd(read_mapped);
            for (i = 0; i < genome_numbers.length; ++i){
                genome = genome_numbers[i];
                num_genomic_mapping[genome].getAndAdd(reads_to_genome[genome]);
            }
            //System.out.println("\nnum_exceptions:"+num_exceptions+" num_neighbors:"+num_neighbors+" num_kmers:"+num_kmers+" num_found:"+num_found);
        }
        
        public void anchor_read(int read_num){
            Node node = null;
            int i, mate, num_segments, step;
            long cur_index, prev_node_id;
            num_segments = paired_end ? 2 : 1;
            for (mate = 0; mate < num_segments; ++mate){
                get_read(read_num, mate);
                anchor_position = reads_scanner[mate].initialize_kmer(0);
                prev_node_id = 0;
                step = reads[mate].length() / KMER_SAMPLES;
                while (anchor_position != Integer.MIN_VALUE && anchor_position <= reads[mate].length() - 31){
                    cur_index = reads_scanner[mate].find_curr_kmer(index_database);
                    try{
                        if (cur_index != -1l){
                            index_database.get_pointer(anchor_pointer, cur_index);
                            if (prev_node_id != anchor_pointer.node_id){ // to save a bit of time
                                node = graphDb.getNodeById(anchor_pointer.node_id);
                                prev_node_id = anchor_pointer.node_id;
                            } 
                            find_locations(node, mate);
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
        
        public void get_read(int read_num, int mate){
            int i;
            reads[mate].clear();
            reads_scanner[mate].set_sequence(read_num);
            reads_scanner[mate].get_complete_sequence(reads[mate].forward_seq, 1, read_num, true);
            reads_scanner[mate].get_complete_sequence(reads[mate].reverse_seq, 1, read_num, false);
            //reads_scanner[mate].get_sequence_quality(reads[mate].quality, 1, read_num);
            reads_scanner[mate].get_sequence_title(reads[mate].name, 1, read_num);
            read_name = reads[mate].name;
            if (paired_end && mate == 1 && reads[0].name.charAt(reads[0].name.length() - 1) == '1')
                read_name.setLength(read_name.length() - 2);
            if (reads[mate].length() > MAX_ALIGNMENT_LENGTH){
                System.out.println(reads[mate].length() + " exceeds the MAX_ALIGNMENT_LENGTH " + MAX_ALIGNMENT_LENGTH +
                        ". Try a larger MAX_ALIGNMENT_LENGTH with --max-length option.");
                System.exit(1);
                
            }
            if ((i = reads[mate].name.indexOf(" ")) != -1)
                reads[mate].name.setLength(i);
            if (DEBUG) 
                System.out.println("read" + read_num);
            /*if (reads[mate].quality.length() > 1){ // input is fastq
                for (; stop >= 0; --stop)
                    if (reads[mate].quality.charAt(stop) > MIN_BASE_QUALITY)
                        break;
                for (; start <= stop; ++start)
                    if (reads[mate].quality.charAt(start) > MIN_BASE_QUALITY)
                        break;
            }*/           
        }
        
        public void find_locations(Node node, int mate) {
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
                    find_address(address, seq_id);
                    genome = address[0];
                    if (locations[mate][genome] != null){
                        sequence = address[1];
                        if (!sequences[mate][genome].contains(sequence))
                            sequences[mate][genome].add(sequence);
                    // calculate the locations based on the offsets in the node    
                        location_array = (int[])r.getProperty(seq_id);
                        seq_len = genome_database.sequence_length[genome][sequence];
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
        
        public void find_address(int[] a, String prp){
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

        public void find_hits(int genome){
            int sequence, num_segments, mate, first_start, prev_start;
            int ref_start = -1, i, n, count;
            num_segments = paired_end ? 2 : 1;
            for (mate = 0; mate < num_segments; ++mate){
                hits[mate].clear();
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
                            if (ref_start - prev_start > MAXBOUND){ 
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
                        //check_add_hit(mate, genome, sequence, hit_counts.get(0)[0]);
                        n = Math.min(hit_counts.size(), 20);//hit_counts.get(0)[1] / 2;
                        for (i = 0; i < n; ++i)
                            check_add_hit(mate, genome, sequence, hit_counts.get(i)[0]);
                        hit_counts.clear();
                        locations[mate][genome][sequence].clear();
                    }
                }
            }
        }
        
        public boolean check_add_hit(int mate, int genome, int sequence, int ref_start){
            int start, stop, offset = 0, bound;
            double score = -1;
            String cigar = null;
            boolean forward = ref_start >= 0;
            hit curr_hit;
            alignment aln;
            if (!forward)
                ref_start = -ref_start - 1;
            for (bound = MINBOUND; bound <= MAXBOUND; bound += MAXBOUND - MINBOUND)
            {
            reference.setLength(0);
            start = ref_start - bound;
            stop = start + reads[mate].length() + 2 * bound - 1;
            //System.out.println(ref_start + " " + genome_database.sequence_length[genome][sequence]);
            if (start >= 0 && stop <= genome_database.sequence_length[genome][sequence] - 1){
                genome_scanner.get_sub_sequence(reference, genome, sequence, start, stop - start + 1, true);
                if ((aln = find_alignment(reference, mate)) != null){
                    score = aln.hit.score;
                    cigar = new String(aln.hit.cigar);
                    offset = aln.hit.offset;
                    curr_hit = new hit(genome, sequence, score, start, offset, forward, cigar);
                    hits[mate].add(curr_hit);
                    return true;
                } else {
                    bounded_aligner.initialize_bound(bound, reads[mate].length());
                    bounded_aligner.align(forward?reads[mate].forward_seq:reads[mate].reverse_seq, reference); 
                    score = bounded_aligner.get_similarity_percentage();
                    if (score >= MIN_ALIGNMENT_SCORE){ // to prevent false positives e.g. random hits
                        cigar = bounded_aligner.get_cigar();
                        offset = bounded_aligner.get_offset();
                        curr_hit = new hit(genome, sequence, score, start, offset, forward, cigar);
                        hits[mate].add(curr_hit);
                        alignments[mate].add(new alignment(reference.toString(),curr_hit));
                        return true;
                    }
                }
            } else { // may map at the borders
                genome_scanner.get_sub_sequence(reference, genome, sequence, ref_start, reads[mate].length(), true);
                aligner.align(forward?reads[mate].forward_seq:reads[mate].reverse_seq, reference); 
                score = aligner.get_similarity_percentage();
                if (score >= MIN_ALIGNMENT_SCORE){ // to prevent false positives e.g. random hits
                    cigar = aligner.get_cigar();
                    //cigar = aligner.get_cigar();
                    offset = aligner.get_offset();
                    curr_hit = new hit(genome, sequence, score, ref_start, offset, forward, cigar);
                    hits[mate].add(curr_hit);
                    alignments[mate].add(new alignment(reference.toString(),curr_hit));
                    return true;
                } 
            }
            }
            return false;
        }

        public alignment find_alignment(StringBuilder reference, int mate) {
            alignment aln = null;
            Iterator<alignment> itr = alignments[mate].iterator();
            while (itr.hasNext()) {
                aln = itr.next();
                if (are_equal(aln.reference, reference))
                    return aln;
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
                
        public boolean single_read_mapped(int genome){
            double first_score;
            hit h;
            boolean secondary;
            if (hits[0].isEmpty()){
                write_single_sam_record(genome, new hit(0, 0, 0, 0, 0, true, null), 0);
                return false;
            } else {
                switch (ALIGNMENT_MODE){
                    case 0:
                        write_single_sam_record(genome, hits[0].remove(), 0);
                        break;
                    case 1:
                        for (secondary = false, first_score = hits[0].peek().score; !hits[0].isEmpty() && first_score - (h = hits[0].remove()).score < 1; secondary = true)
                            write_single_sam_record(genome, h, secondary?256:0);
                        break;
                }
                return true;
            }
        }
        
        public boolean paired_read_mapped(int genome){
            double score;
            hit h1, h2, null_hit, best_hit1 = null, best_hit2 = null;
            int best_cigar_length = 0, cigar_length;
            Iterator<hit> itr1, itr2;
            double best_score = -1;
            boolean first_hit = true, had_hit;
            best_score = 0;
            if (hits[0].isEmpty() || hits[1].isEmpty()){
                null_hit = new hit(0, 0, 0, 0, 0, true, null);
                if (hits[0].isEmpty() && hits[1].isEmpty())
                // multiple segments|each segment properly aligned|secondary|unmapped
                    write_paired_sam_record(genome, null_hit, null_hit, 1|0|0|12, 1|0|0|12);
                else if (!hits[1].isEmpty())
                    write_paired_sam_record(genome, null_hit, hits[1].remove(), 1|0|0|4, 1|0|0|8);
                else 
                    write_paired_sam_record(genome, hits[0].remove(), null_hit, 1|0|0|8, 1|0|0|4);
                return false;
            } else {
                switch (ALIGNMENT_MODE){
                    case 0:
                        for (had_hit = false, itr1 = hits[0].iterator(); itr1.hasNext(); ){
                            h1 = itr1.next();
                            for (itr2 = hits[1].iterator(); itr2.hasNext(); ){
                                h2 = itr2.next();
                                //System.out.println(h1.toString() + " & " +h2.toString());
                                if (h1.sequence == h2.sequence && h1.forward == !h2.forward && properly_mapped(h1, h2)){
                                    if (first_hit){
                                        first_hit = false;
                                        best_hit1 = h1;
                                        best_hit2 = h2;
                                        best_score = (h1.score + h2.score) / 2;
                                        best_cigar_length = h1.cigar.length() + h1.cigar.length();
                                        had_hit = true;
                                    } else{ 
                                        score = (h1.score + h2.score) / 2;
                                        cigar_length = h1.cigar.length() + h1.cigar.length();
                                        if (best_score < score || 
                                           (best_score == score && cigar_length  < best_cigar_length) || 
                                           (best_score == score && cigar_length == best_cigar_length && h1.start < best_hit1.start)){
                                            best_hit1 = h1;
                                            best_hit2 = h2;
                                            best_score = score;
                                            best_cigar_length = cigar_length;
                                            had_hit = true;
                                        }
                                    }
                                }
                                
                            }
                        } 
                        if (!had_hit)
                            write_paired_sam_record(genome, hits[0].remove(), hits[1].remove(), 1|0|0|0, 1|0|0|0);
                        else
                            write_paired_sam_record(genome, best_hit1, best_hit2, 1|2|0|0, 1|2|0|0);
                        break;
                    case 1:
                    // finds best score    
                        for (had_hit = false, itr1 = hits[0].iterator(); itr1.hasNext(); ){
                            h1 = itr1.next();
                            for (itr2 = hits[1].iterator(); itr2.hasNext(); ){
                                h2 = itr2.next();
                                if (h1.sequence == h2.sequence && h1.forward == !h2.forward && properly_mapped(h1, h2)){
                                    if (first_hit){
                                        first_hit = false;
                                        best_score = (h1.score + h2.score) / 2;
                                        had_hit = true;
                                    } else{ 
                                        score = (h1.score + h2.score) / 2;
                                        if (best_score - score < 1){
                                            had_hit = true;
                                        }
                                    }
                                }
                                
                            }
                        } 
                        if (!had_hit)
                            write_paired_sam_record(genome, hits[0].remove(), hits[1].remove(), 1|0|0|0, 1|0|0|0);
                        else{ // writes best_scored hits
                            for (itr1 = hits[0].iterator(); itr1.hasNext(); ){
                                h1 = itr1.next();
                                for (itr2 = hits[1].iterator(); itr2.hasNext(); ){
                                    h2 = itr2.next();
                                    if (h1.sequence == h2.sequence && h1.forward == !h2.forward && properly_mapped(h1, h2)){
                                        if (first_hit){
                                            score = (h1.score + h2.score) / 2;
                                            if (best_score - score < 1){
                                                first_hit = false;
                                                write_paired_sam_record(genome, h1, h2, 1|2|0|0, 1|2|0|0);
                                            }
                                        } else{ 
                                            score = (h1.score + h2.score) / 2;
                                            if (best_score - score < 1){
                                                write_paired_sam_record(genome, h1, h2, 1|2|256|0, 1|2|256|0);
                                            }
                                        }
                                    }

                                }
                            } 
                        }
                        break;
                }
                return true;
            }    
        }
        
        public boolean properly_mapped(hit h1, hit h2){
            int position1 = h1.start + h1.offset;
            int position2 = h2.start + h2.offset;
            if (h1.forward)
                return position2 - (position1 + reads[0].forward_seq.length() ) >= -reads[0].forward_seq.length() && position2 - (position1 + reads[0].forward_seq.length() ) < (2 * INNER_READS_DISTANCE) ;
            else
                return position1 - (position2 + reads[1].forward_seq.length() ) >= -reads[1].forward_seq.length() && position1 - (position2 + reads[1].forward_seq.length() ) < (2 * INNER_READS_DISTANCE);
        }

        public void write_single_sam_record(int genome, hit h, int secondary){
            int flag = secondary;
            sam_record.setLength(0);
            if (h.genome == 0){
                flag = 4; //unmapped
                sam_record.append(read_name).append('\t');
                sam_record.append(flag).append('\t');
                sam_record.append('*').append('\t');
                sam_record.append(0).append('\t');
                sam_record.append(0).append('\t');
                sam_record.append('*').append('\t');
                sam_record.append('*').append('\t');
                sam_record.append('0').append('\t');
                sam_record.append('0').append('\t');
                sam_record.append(reads[0].forward_seq).append('\t');
                sam_record.append('*').append('\n');
                //sam_record.append(reads[0].quality).append('\n');
                
            } else {
                flag = h.forward?0:16; // direction
                sam_record.append(reads[0].name).append('\t');
                sam_record.append(flag).append('\t');
                sam_record.append(genome_database.sequence_titles[h.genome][h.sequence].split("\\s")[0]).append('\t');
                sam_record.append(h.start + h.offset + 1).append('\t');
                sam_record.append((int)Math.round(-10 * Math.log10(1 - (h.score - 1) / 100.0))).append('\t');
                sam_record.append(h.cigar).append('\t');
                sam_record.append('*').append('\t');
                sam_record.append('0').append('\t');
                sam_record.append('0').append('\t');
                sam_record.append(h.forward?reads[0].forward_seq:reads[0].reverse_seq).append('\t');
                sam_record.append('*').append('\n');
                //sam_record.append(reads[0].quality).append('\n');
            }
            try {
                mapped[genome].write(sam_record.toString());
            } catch (IOException ex) {
                System.err.println(ex.getMessage());
            }
        }

        public void write_paired_sam_record(int genome, hit h1, hit h2, int flag1, int flag2){
        // the first segment    
            flag1 |= 64; 
            sam_record.setLength(0);
            int position1 = h1.start + h1.offset + 1;
            int position2 = h2.start + h2.offset + 1;
            /*if (position1 < 1)
                correct_initial_hit(h1);
            if (position2 > genomeDb.sequence_length[h2.genome][h2.sequence])
                correct_final_hit(h2);*/
            if (h1.genome == 0){
                sam_record.append(read_name).append('\t');
                sam_record.append(flag1).append('\t');
                sam_record.append('*').append('\t');
                sam_record.append(0).append('\t');
                sam_record.append(0).append('\t');
                sam_record.append('*').append('\t');
                sam_record.append('*').append('\t');
                sam_record.append('0').append('\t');
                sam_record.append('0').append('\t');
                sam_record.append(reads[0].forward_seq).append('\t');
                sam_record.append('*').append('\n');
                //sam_record.append(reads[0].quality).append('\n');
                
            } else {
                flag1 |= h1.forward?0:16; // SEQ being reverse complemented
                flag1 |= h1.forward?32:0; // SEQ being reverse complemented
            //qname
                sam_record.append(read_name).append('\t'); 
            //flag    
                sam_record.append(flag1).append('\t'); 
            //rname    
                sam_record.append(genome_database.sequence_titles[h1.genome][h1.sequence].split("\\s")[0]).append('\t');
            //pos    
                sam_record.append(position1).append('\t');
            //mapq    
                sam_record.append((int)Math.round(-10 * Math.log10(1 - (h1.score/100.0 - 0.01)))).append('\t');
            //cigar    
                sam_record.append(h1.cigar).append('\t');
            //rnext
                if (h1.sequence == h2.sequence)
                    sam_record.append('=').append('\t');
                else if (h2.genome == 0)
                    sam_record.append('*').append('\t');
                else 
                    sam_record.append(genome_database.sequence_titles[h2.genome][h2.sequence].split("\\s")[0]).append('\t');
            //pnext    
                sam_record.append(position2).append('\t');
            //tlen
                if (position2 == 0)
                    sam_record.append('0').append('\t');
                else if (position2 > position1)
                    sam_record.append(position2 - position1 + reads[0].forward_seq.length()).append('\t');
                else
                    sam_record.append('-').append(position1 - position2 + reads[1].forward_seq.length()).append('\t');
            //seq    
                sam_record.append(h1.forward?reads[0].forward_seq:reads[0].reverse_seq).append('\t');
            //qual    
                sam_record.append('*').append('\n');
                //sam_record.append(reads[0].quality).append('\n');
            }
        // the second segment    
            flag2 |= 128; 
            if (h2.genome == 0){ // not mapped
                sam_record.append(read_name).append('\t');
                sam_record.append(flag2).append('\t');
                sam_record.append('*').append('\t');
                sam_record.append(0).append('\t');
                sam_record.append(0).append('\t');
                sam_record.append('*').append('\t');
                sam_record.append('*').append('\t');
                sam_record.append('0').append('\t');
                sam_record.append('0').append('\t');
                sam_record.append(reads[1].forward_seq).append('\t');
                sam_record.append('*').append('\n');
                //sam_record.append(reads[1].quality).append('\n');
                
            } else {
                flag2 |= h2.forward?0:16; // SEQ of the next segment in the template being reverse complemented
                flag2 |= h2.forward?32:0; // SEQ of the next segment in the template being reverse complemented
                sam_record.append(read_name).append('\t');
                sam_record.append(flag2).append('\t');
                sam_record.append(genome_database.sequence_titles[h2.genome][h2.sequence].split("\\s")[0]).append('\t');
                sam_record.append(position2).append('\t');
                sam_record.append((int)Math.round(-10 * Math.log10(1 - (h2.score/100.0 - 0.01)))).append('\t');
                sam_record.append(h2.cigar).append('\t');
                if (h1.sequence == h2.sequence)
                    sam_record.append('=').append('\t');
                else if (h1.genome == 0)
                    sam_record.append('*').append('\t');
                else 
                    sam_record.append(genome_database.sequence_titles[h1.genome][h1.sequence].split("\\s")[0]).append('\t');
                sam_record.append(position1).append('\t');
                if (position1 == 0)
                    sam_record.append('0').append('\t');
                else if (position2 > position1)
                    sam_record.append('-').append(position2 - position1 + reads[0].forward_seq.length()).append('\t');
                else
                    sam_record.append(position1 - position2 + reads[1].forward_seq.length()).append('\t');
                sam_record.append(h2.forward?reads[1].forward_seq:reads[1].reverse_seq).append('\t');
                sam_record.append('*').append('\n');
                //sam_record.append(reads[1].quality).append('\n');
            }
            try {
               mapped[genome].write(sam_record.toString());
            } catch (IOException ex) {
                System.err.println(ex.getMessage());
            }
        }
        
        /*void correct_initial_hit(hit h){
            int i, num = 0, len, pos = h.start + h.offset + 1;
            char ch, last_operation;
            StringBuilder cigar = new StringBuilder();
            len = h.cigar.length();
            for (i = 0; i < len && pos < 1; ++i){
                ch = h.cigar.charAt(i);
                if (ch >= 48 && ch < 58)
                    num = num * 10 + ch;
                else {
                    if (ch == 'S')
                        h.
                    if (pos + num <)
                    pos += num;
                    h
                    num = 0;
                }
                    
            }
        }*/

    }   
    
    /**
     * Constructs a pangenome database from given genomes.
     * 
     * @param genome_paths_file Path to the FASTA genome files. 
     * @param pangenome_path Path to the database folder
     */  
    public void initialize_pangenome() {
    // If a database folder is already exist in the specified path, removes all the content of it.
        Node pangenome_node;
        if (PATH_TO_THE_PANGENOME_DATABASE == null){
            System.out.println("PATH_TO_THE_PANGENOME_DATABASE is empty.");
            System.exit(1);
        }
        File theDir = new File(PATH_TO_THE_PANGENOME_DATABASE);
        if (theDir.exists()) {
            try {
                FileUtils.deleteRecursively(new File(PATH_TO_THE_PANGENOME_DATABASE));
            } catch (IOException ioe) {
                System.out.println("Failed to delete the database " + PATH_TO_THE_PANGENOME_DATABASE);
                System.exit(1);  
            }
        } else {
            try {
                theDir.mkdir();
            } catch (SecurityException se) {
                System.out.println("Failed to create directory " + PATH_TO_THE_PANGENOME_DATABASE);
                System.exit(1);
            }
        }
        if (PATH_TO_THE_GENOMES_FILE == null){
            System.out.println("PATH_TO_THE_GENOMES_FILE is empty.");
            System.exit(1);
        }
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();  
        registerShutdownHook(graphDb);
        startTime = System.currentTimeMillis();
        num_nodes = 0;
        num_edges = 0;
        num_bases = 0;
        num_degenerates = 0;
        genomeDb = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH, PATH_TO_THE_GENOMES_FILE);
        indexDb = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH, PATH_TO_THE_GENOMES_FILE, genomeDb, K_SIZE);
        K_SIZE = indexDb.get_K();
        genome_scanner = new SequenceScanner(genomeDb, 1, 1, K_SIZE, indexDb.get_pre_len());
        try (Transaction tx = graphDb.beginTx()) {
            pangenome_node = graphDb.createNode(pangenome_label);
            pangenome_node.setProperty("k_mer_size", K_SIZE);
            pangenome_node.setProperty("date", new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(new Date()));
            tx.success();
        }
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
        int j, len, previous_num_genomes;
        long byte_number = 0;
        int[] address = new int[4];
        Node pangenome_node, start, seq_node;
        if (PATH_TO_THE_GENOMES_FILE == null){
            System.out.println("PATH_TO_THE_GENOMES_FILE is empty.");
            System.exit(1);
        }
        if (! new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH).exists()) {
            System.out.println("No database found in " + PATH_TO_THE_PANGENOME_DATABASE);
            System.exit(1);
        }
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
        // if the genome database is not available, reconstruct it.    
            if (!Files.exists(Paths.get(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH))) {
            // read genomes information from the graph and rebuild the genomes database
                genomeDb = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH, graphDb);
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
                            genomeDb.genomes_buff[(int) (byte_number / genomeDb.parts_size[0])].put((byte) ((genomeDb.binary[seq.charAt(j)] << 4) | genomeDb.binary[seq.charAt(j + 1)]));
                        }
                        if (len == seq.length() - 1) {
                            genomeDb.genomes_buff[(int) (byte_number / genomeDb.parts_size[0])].put((byte) (genomeDb.binary[seq.charAt(len)] << 4));
                            ++byte_number;
                        }
                    }
                }
            } else {
                genomeDb = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH);
            }
            genomeDb.add_genomes(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH, PATH_TO_THE_GENOMES_FILE);
            indexDb = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH, PATH_TO_THE_GENOMES_FILE, genomeDb, graphDb, previous_num_genomes);
            genome_scanner = new SequenceScanner(genomeDb, previous_num_genomes + 1, 1, K_SIZE, indexDb.get_pre_len());
            tx.success();
        }
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

    public void remove_genomes() {
        
    }
    
    public void map_reads() {
        int i, j, genome, n = 0;
        Node pangenome_node;
        int[] genome_numbers = null;
        String line;
        BufferedReader in;
        boolean paired;
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
        sequencingDb = new SequenceDatabase[2];
        sequencingDb[0] = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + READS_DATABASE_PATH + ".1", PATH_TO_THE_FIRST_SRA);
        if (paired){
            sequencingDb[1] = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + READS_DATABASE_PATH + ".2", PATH_TO_THE_SECOND_SRA);
            if (sequencingDb[0].num_sequences[1] != sequencingDb[1].num_sequences[1]){
                System.err.println("Paired-end libraries contain different number of reads.");
                System.exit(1);
            }
        }
        genomeDb = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH);
        num_genomic_mapping = new AtomicInteger[genomeDb.num_genomes + 1];
        num_mapping = new AtomicInteger(0);
        read_number = new AtomicInteger(1);
        //number_of_alignments = new AtomicInteger(0);
        try {
            in = new BufferedReader(new FileReader(PATH_TO_THE_GENOME_NUMBERS_FILE));
            while (in.ready())
                if (!in.readLine().trim().equals(""))
                    ++n;
            in.close();
            genome_numbers = new int[n];
            in = new BufferedReader(new FileReader(PATH_TO_THE_GENOME_NUMBERS_FILE));
            for (i = 0; in.ready(); ){
                line = in.readLine().trim();
                if (!line.equals(""))
                    genome_numbers[i++] = Integer.parseInt(line);
            }
            Arrays.sort(genome_numbers);
            in.close();
        } catch (Exception ex){
            System.err.println("Error in reading genome numbers");
        }
        System.out.println("Mapping " + sequencingDb[0].num_sequences[1] + " reads on " + genome_numbers.length + " genome(s) :");
        System.out.print("0..................................................100\n ");
        BufferedWriter[] out = new BufferedWriter[genomeDb.num_genomes + 1];
        try{
            for (i = 0; i < genome_numbers.length; ++i){
                genome = genome_numbers[i];
                out[genome] = new BufferedWriter(new FileWriter(
                        OUTPUT_PATH + "/pantools_"+ genome + ".sam"));
                out[genome].write("@HD\tVN:1.5\tSO:unsorted\n");
                for (j = 1; j <= genomeDb.num_sequences[genome]; ++j){
                    out[genome].write("@SQ\tSN:" + genomeDb.sequence_titles[genome][j].split("\\s")[0]  + "\tLN:" + genomeDb.sequence_length[genome][j] + "\n");
                }
                out[genome].write("@PG\tID:PanTools\tPN:PanTools\n");
                num_genomic_mapping[genome] = new AtomicInteger(0);
            }
            genomeDb.close();
            try{
                ExecutorService es = Executors.newFixedThreadPool(THREADS);
                for (i = 0; i < THREADS; ++i)
                    es.execute(new Map(i, genome_numbers, out, paired));
                es.shutdown();
                es.awaitTermination(10, TimeUnit.DAYS);        
            } catch (InterruptedException e){

            }
            new File(PATH_TO_THE_PANGENOME_DATABASE + READS_DATABASE_PATH + ".1/sequences.db").delete();
            new File(PATH_TO_THE_PANGENOME_DATABASE + READS_DATABASE_PATH + ".2/sequences.db").delete();
            new File(PATH_TO_THE_PANGENOME_DATABASE + READS_DATABASE_PATH + ".1/sequences.info").delete();
            new File(PATH_TO_THE_PANGENOME_DATABASE + READS_DATABASE_PATH + ".2/sequences.info").delete();
            new File(PATH_TO_THE_PANGENOME_DATABASE + READS_DATABASE_PATH + ".1").delete();
            new File(PATH_TO_THE_PANGENOME_DATABASE + READS_DATABASE_PATH + ".2").delete();
            System.out.println("\nTotal mapping = " + num_mapping);
            for (i = 0; i < genome_numbers.length; ++i){
                genome = genome_numbers[i];
                out[genome].close();
                System.out.println(num_genomic_mapping[genome] + " mapped to genome " + genome);
            }
            //System.out.println("Number of pairwise alignments : " + number_of_alignments);
        } catch(IOException ex){

        }
        graphDb.shutdown();
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
        indexDb = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH);
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
                        start_ptr = locate(address);
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
        indexDb = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH);
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
                    start = locate(new int[]{genome1, seq1, loc1 + 1});
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
    public IndexPointer locate(int[] addr) {
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
                System.out.println("Processing genome " + genome_scanner.get_genome() + " :             ");
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
