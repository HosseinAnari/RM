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
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import sequence.SequenceDatabase;
import sequence.SequenceScanner;
import index.IndexPointer;
import index.IndexDatabase;
import index.IndexScanner;
import index.kmer;
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
import java.util.HashSet;
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
import static pantools.Pantools.nucleotide_label;
import static pantools.Pantools.pangenome_label;
import static pantools.Pantools.sequence_label;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.neo4j.graphdb.factory.GraphDatabaseSettings;
import static pantools.Pantools.ALIGNMENT_BOUND;
import static pantools.Pantools.ALIGNMENT_MODE;
import static pantools.Pantools.INTERLEAVED;
import static pantools.Pantools.NUM_KMER_SAMPLES;
import static pantools.Pantools.MAX_ALIGNMENT_LENGTH;
import static pantools.Pantools.MAX_FRAGMENT_LENGTH;
import static pantools.Pantools.MIN_HIT_LENGTH;
import static pantools.Pantools.MIN_IDENTITY;
import static pantools.Pantools.MAX_NUM_LOCATIONS;
import static pantools.Pantools.SHOULDER;
import static pantools.Pantools.binary;
import static pantools.Pantools.genomeDb;
import static pantools.Pantools.genomeSc;
import static pantools.Pantools.get_lines_count;
import static pantools.Pantools.graphDb;
import static pantools.Pantools.heapSize;
import static pantools.Pantools.indexDb;
import static pantools.Pantools.indexSc;
import static pantools.Pantools.open_file;
import static pantools.Pantools.reverse_complement;

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
    private SequenceScanner genomeSc;    
    private BlockingQueue<FastqRecord>[][] fastq_records; 
    
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
        public int score;
        public int start;
        public int offset;
        public int length;
        public int deletions;
        public boolean forward;
        public String cigar;
        public String reference;
        public single_hit(int gn, int sq, double idn, int sc, int st, int off, int len, int del, boolean fw, String cg, String r){
            genome = gn;
            sequence = sq;
            identity = idn;
            score = sc;
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
            score = h.score;
            start = h.start;
            offset = h.offset;
            length = h.length;
            deletions = h.deletions;
            forward = h.forward;
            cigar = h.cigar;
            reference = h.reference;
        }
        
        public single_hit(){
        }
        

        @Override
        public String toString(){
            return "(genome:" + genome + 
                    ",sequence:" + sequence + 
                    ",identity:" + identity + 
                    ",score:" + score + 
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
        public single_hit h1;
        public single_hit h2;
        public paired_hit(int flen, int gn1, int sq1, double idn1, int sc1, int st1, int off1, int len1, int del1, boolean fw1, String cg1, String r1,
                          int gn2, int sq2, double idn2, int sc2, int st2, int off2, int len2, int del2, boolean fw2, String cg2, String r2){
            fragment_length = flen;
            h1.genome = gn1;
            h1.sequence = sq1;
            h1.identity = idn1;
            h1.score = sc1;
            h1.start = st1;
            h1.offset = off1;
            h1.length = len1;
            h1.deletions = del1;
            h1.forward = fw1;
            h1.cigar = cg1;
            h1.reference = r1;

            h2.genome = gn2;
            h2.sequence = sq2;
            h2.identity = idn2;
            h2.score = sc2;
            h2.start = st2;
            h2.offset = off2;
            h2.length = len2;
            h2.deletions = del2;
            h2.forward = fw2;
            h2.cigar = cg2;
            h2.reference = r2;
        }

        public paired_hit(int f_len, single_hit hit1, single_hit hit2){
            fragment_length = f_len;
            h1 = hit1;
            h2 = hit2;
        }

        public int get_score(){
            return h1.score + h2.score;
        }

        public int get_min_start(){
            return Math.min(h1.start, h2.start);
        }

        public int get_max_start(){
            return Math.max(h1.start, h2.start);
        }

        @Override
        public String toString(){
            return "(genome1:" + h1.genome + 
                    ",sequence1:" + h1.sequence + 
                    ",identity1:" + h1.identity + 
                    ",score1:" + h1.score + 
                    ",start1:" + h1.start + 
                    ",offset1:" + h1.offset + 
                    ",length1:" + h1.length + 
                    ",deletions1:" + h1.deletions + 
                    ",forward1:" + h1.forward + 
                    ",reference1:" + h1.reference + 
                    ",cigar1:" + h1.cigar +")" +
                    "\n" +
                    "(genome2:" + h2.genome + 
                    ",sequence2:" + h2.sequence + 
                    ",identity2:" + h2.identity + 
                    ",score2:" + h2.score + 
                    ",start2:" + h2.start + 
                    ",offset2:" + h2.offset + 
                    ",length2:" + h2.length + 
                    ",deletions2:" + h2.deletions + 
                    ",forward2:" + h2.forward + 
                    ",reference2:" + h2.reference + 
                    ",cigar2:" + h2.cigar +")";
        }
    }

    /**
     * Implements a comparator for integer arrays of size two
     */
    public static class single_hitComparator implements Comparator<single_hit> {
        @Override
        public int compare(single_hit x, single_hit y) {
            if (x.score > y.score) 
                return -1;
            else if (x.score < y.score) 
                return 1;
            else if (x.sequence > y.sequence) 
                return -1;
            else if (x.sequence < y.sequence) 
                return 1;
            else if (x.start > y.start) 
                return -1;
            else if (x.start < y.start) 
                return 1;
            else 
                return 0;
        }
    }      

    public static class paired_hitComparator implements Comparator<paired_hit> {
        @Override
        public int compare(paired_hit x, paired_hit y) {
            if (x.get_score() > y.get_score()) 
                return -1;
            else if (x.get_score() < y.get_score()) 
                return 1;
            else if (x.fragment_length < y.fragment_length) 
                return -1;
            else if (x.fragment_length > y.fragment_length) 
                return 1;
            else 
                return 0;
        }
    }      

    public static class IntPairComparator implements Comparator<int[]> {
        @Override
        public int compare(int[] x, int[] y) {
            if (x[0] > y[0]) 
                return -1;
            else if (x[0] < y[0]) 
                return 1;
            else if (x[1] < y[1]) 
                return -1;
            else if (x[1] > y[1]) 
                return 1;
            else
                return 0;
        }
    }    

    public static class IntComparator implements Comparator<Integer> {

        @Override
        public int compare(Integer v1, Integer v2) {
            return v1 < v2 ? -1 : v1 > v2 ? 1 : 0;
        }
    }

    public class Generate_fatsq_records implements Runnable {
        boolean paired;
        public Generate_fatsq_records(boolean p){
            paired = p;
        }
        @Override
        public void run() {
            int t;
            long counter;
            FastqReader[] reader = new FastqReader[2];
            reader[0] = new FastqReader(new File(PATH_TO_THE_FIRST_SRA), true);
            if (paired && !INTERLEAVED)
                reader[1] = new FastqReader(new File(PATH_TO_THE_SECOND_SRA), true);
            try {
                for (counter = 0; reader[0].hasNext(); ++counter){
                    t = (int)(counter % THREADS);
                    fastq_records[t][0].put(reader[0].next());
                    if (paired){ 
                        if(INTERLEAVED)
                            fastq_records[t][1].put(reader[0].next());
                        else
                            fastq_records[t][1].put(reader[1].next());
                    }
                    if (counter % 1000 == 0)
                        System.out.print("\rReads: " + counter);
                }
                System.out.println("\rReads: " + counter);
                for (t = 0; t < THREADS; ++t){
                   fastq_records[t][0].put(new FastqRecord("", "", "", "")); 
                   if (paired)
                       fastq_records[t][1].put(new FastqRecord("", "", "", "")); 
                }
            } catch (InterruptedException ex) {
                System.err.println(ex.getMessage());
                Logger.getLogger(GenomeLayer.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    } 
    
    public class Map implements Runnable {
        IndexScanner indexSc;
        SequenceScanner genomeSc;
        int thread_id;
        int K;
        FastqRecord[] reads;
        ArrayList<int[]>[][]locations;
        IntComparator intcomp = new IntComparator();
        IndexPointer pointer;
        StringBuilder reference;
        BoundedLocalSequenceAlignment bounded_aligner;
        LocalSequenceAlignment aligner;
        PriorityQueue<single_hit>[][] hits;
        LinkedList<single_hit>[] alignments;
        PriorityQueue<single_hit> single_hits;
        Queue<single_hit> single_hits_2;
        PriorityQueue<paired_hit> paired_hits;
        Queue<paired_hit> paired_hits_2;
        boolean paired_end;
        IntPairComparator int_pair_comp;
        ArrayList<int[]> hit_counts;
        int num_neighbors = 0;
        int num_kmers = 0;
        int num_found = 0;
        int[] address = new int[3];
        ArrayList<Integer> genome_numbers;
        SAMFileWriter[] sam_writers;
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
        paired_hitComparator phc = new paired_hitComparator();
        single_hitComparator shc = new single_hitComparator();
        StringBuilder[] forward_read;
        StringBuilder[] reverse_read;
        StringBuilder[] quality;
        kmer current_kmer; 
        int[] read_len;
        String[] read_name;
        int num_hits = 0;
        int num_alns = 0;
        single_hit alignment_result;
        HashSet<Integer>[] unmapped_genomes;
        single_hit[] best_hit;
        int num_segments;
        
        public Map(int id, ArrayList<Integer> gn, boolean paired, SAMFileWriter[] sams) {
            int i, j, genome, abun;
            indexSc = new IndexScanner(indexDb);
            genomeSc = new SequenceScanner(genomeDb, 1, 1, indexSc.get_K(), indexSc.get_pre_len());
            K = indexSc.get_K();
            current_kmer = new kmer(indexSc.get_K(), indexSc.get_pre_len());
            num_genomes = genomeDb.num_genomes;
            unmapped_genomes = new HashSet[2];
            num_sequences = new int[num_genomes + 1];
            sequence_length = new long[num_genomes + 1][];
            sequence_titles = new String[num_genomes + 1][];
            sam_writers = sams;
            thread_id = id;
            genome_numbers = gn;
            genome_numbers.sort(intcomp);
            paired_end = paired; 
            num_segments = paired_end ? 2 : 1;
            alignment_result = new single_hit();
            hits = new PriorityQueue[2][];
            alignments = new LinkedList[2];
            locations = new ArrayList[2][];
            reads = new FastqRecord[2];
            quality = new StringBuilder[2];
            forward_read = new StringBuilder[2];
            reverse_read = new StringBuilder[2];
            read_len = new int[2];
            read_name = new String[2];
            best_hit = new single_hit[2];
            
            quality[0] = new StringBuilder();
            forward_read[0] = new StringBuilder();
            reverse_read[0] = new StringBuilder();
            hits[0] = new PriorityQueue[num_genomes + 1];
            alignments[0] = new LinkedList();
            single_hits = new PriorityQueue(shc);
            single_hits_2 = new LinkedList();
            paired_hits = new PriorityQueue(phc);
            paired_hits_2 = new LinkedList();
            locations[0] = new ArrayList[num_genomes + 1];
            unmapped_genomes[0] = new HashSet();
            genome_sizes = new long[num_genomes + 1];
            shared = new int[num_genomes + 1];
            unique = new int[num_genomes + 1];
            unmapped = new int[num_genomes + 1];
            raw_abundance = new double[num_genomes + 1];
            for (i = 0; i <= num_genomes; ++ i){
                hits[0][i] = null;
                locations[0][i] = null;
                genome_sizes[i] = genomeDb.genome_length[i];
                num_sequences[i] = genomeDb.num_sequences[i];
                sequence_length[i] = new long[num_sequences[i] + 1];
                sequence_titles[i] = new String[num_sequences[i] + 1];
                shared[i] = 0;
            }
            for (i = 0; i < genome_numbers.size(); ++i){
                genome = genome_numbers.get(i);
                hits[0][genome] = new PriorityQueue(shc);
                locations[0][genome] = new ArrayList();
                for (j = 1; j <= num_sequences[genome]; ++j){
                    sequence_length[genome][j] = genomeDb.sequence_length[genome][j];
                    sequence_titles[genome][j] = genomeDb.sequence_titles[genome][j];
                }
            }
            if (paired_end){
                quality[1] = new StringBuilder();
                forward_read[1] = new StringBuilder();
                reverse_read[1] = new StringBuilder();
                hits[1] = new PriorityQueue[num_genomes + 1];
                alignments[1] = new LinkedList();
                locations[1] = new ArrayList[num_genomes + 1];
                unmapped_genomes[1] = new HashSet();
                for (i = 0; i <= num_genomes; ++ i){
                    hits[1][i] = null;
                    locations[1][i] = null;
                }
                for (i = 0; i < genome_numbers.size(); ++i){
                    genome = genome_numbers.get(i);
                    hits[1][genome] = new PriorityQueue(shc);
                    locations[1][genome] = new ArrayList();
                }
            }
            pointer = new IndexPointer();
            reference = new StringBuilder();
            //bounded_aligner = new LocalSequenceAlignment(GAP_OPEN, GAP_EXT, MAX_ALIGNMENT_LENGTH, CLIP, 'N');
            bounded_aligner = new BoundedLocalSequenceAlignment(GAP_OPEN, GAP_EXT, MAX_ALIGNMENT_LENGTH, ALIGNMENT_BOUND, CLIPPING_STRINGENCY, 'N');
            aligner = new LocalSequenceAlignment(GAP_OPEN, GAP_EXT, MAX_ALIGNMENT_LENGTH, CLIPPING_STRINGENCY, 'N');
            int_pair_comp = new IntPairComparator();
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
                    System.out.println(ex.getMessage());
                }
            }
        }

        @Override
        public void run() {
            int i, mate, genome;
            Iterator<Integer> itr;
            try{
                reads[0] = fastq_records[thread_id][0].take();
                if (paired_end)
                    reads[1] = fastq_records[thread_id][1].take();
                try(Transaction tx = graphDb.beginTx()){
                    while (reads[0].getReadLength() != 0){
                        get_read();
                        for (mate = 0; mate < num_segments; ++mate){
                            find_locations(mate);
                            itr = genome_numbers.iterator()
;                            while (itr.hasNext()){
                                genome = itr.next();
                                cluster_and_examine_locations(mate, genome, 0);
                            }
                        }
                        report_hits();
                        reads[0] = fastq_records[thread_id][0].take();
                        if (paired_end)
                            reads[1] = fastq_records[thread_id][1].take();
                    }
                    tx.success();
                }
            }catch(InterruptedException e){
                System.err.println(e.getMessage());
            }
            number_of_hits.getAndAdd(num_hits);
            number_of_alignments.getAndAdd(num_alns);
            for (i = 0; i < genome_numbers.size(); ++i){
                genome = genome_numbers.get(i);
                num_shared_mapping[genome].getAndAdd(shared[genome]);
                num_unique_mapping[genome].getAndAdd(unique[genome]);
                num_unmapped[genome].getAndAdd(unmapped[genome]);
            }
        }
        
        public void find_locations(int mate){
            Node node;
            int i, step, pos;
            long cur_index, prev_node_id;
            prev_node_id = 0;
            step = read_len[mate] / NUM_KMER_SAMPLES;
            step = (step == 0 ? 1 : step);
            initialize_kmer(forward_read[mate]);
            for (pos = K; pos < read_len[mate];){
                cur_index = indexSc.find(current_kmer);
                try{
                    if (cur_index != -1l){
                        indexSc.get_pointer(pointer, cur_index);
                        if (prev_node_id != pointer.node_id){ // to save a bit of time
                            node = graphDb.getNodeById(pointer.node_id);
                            prev_node_id = pointer.node_id;
                            explore_node(node, mate, pos - K);
                        } 
                    }
                    if (pos + step >= read_len[mate])
                        break;
                    for (i = 0; i < step; ++i, ++pos)
                        current_kmer.next_kmer(binary[forward_read[mate].charAt(pos)] & 3);
                } catch (NotFoundException|ClassCastException ex){
                    //num_exceptions++;
                    //System.out.println(ex.getMessage());
                } 
                //System.out.println(current_kmer.toString());
            }
        }
        
        public void get_read(){
            read_len[0] = reads[0].getReadLength();
            quality[0].setLength(0);
            quality[0].append(reads[0].getBaseQualityString());
            forward_read[0].setLength(0);
            forward_read[0].append(reads[0].getReadString());
            reverse_read[0].setLength(0);
            reverse_read[0].append(reads[0].getReadString());
            reverse_complement(reverse_read[0]);
            read_name[0] = getBaseId(reads[0].getReadName());
            if (paired_end){
                quality[1].setLength(0);
                quality[1].append(reads[1].getBaseQualityString());
                read_len[1] = reads[1].getReadLength();
                forward_read[1].setLength(0);
                forward_read[1].append(reads[1].getReadString());
                reverse_read[1].setLength(0);
                reverse_read[1].append(reads[1].getReadString());
                reverse_complement(reverse_read[1]);
                read_name[1] = getBaseId(reads[1].getReadName());
            }
        }

        public  String getBaseId(String Id) {
            int slashIdx = Id.indexOf("/");
            int spaceIdx = Id.indexOf(" ");

            if ((slashIdx == -1) && (spaceIdx == -1)) {
                return Id;
            }

            int idx = -1;
            if (slashIdx == -1) {
                idx = spaceIdx;
            } else if (spaceIdx == -1) {
                idx = slashIdx;
            } else {
                idx = spaceIdx < slashIdx ? spaceIdx : slashIdx;
            }

            return Id.substring(0, idx);
        }        
        
        public void initialize_kmer(StringBuilder read){
            int i;
            current_kmer.reset();
            for (i = 0; i < K; ++i)
                current_kmer.next_kmer(binary[read.charAt(i)] & 3);
        }

        public void explore_node(Node node, int mate, int position) {
            int i, loc, offset;
            long seq_len;
            char side;
            int[] location_array;
            int genome, sequence;
            boolean is_canonical;
            long[] frequenceis;
            offset = pointer.offset;
            is_canonical = current_kmer.get_canonical();
            int node_len;
            frequenceis = (long[])node.getProperty("frequencies");
        // for each incoming edge to the node of the anchor    
            for (Relationship r: node.getRelationships(Direction.INCOMING)){
                //num_neighbors++;
                side = r.getType().name().charAt(1);
            // for all seuences passing that node     
                for (String seq_id: r.getPropertyKeys()){
                    extract_address(address, seq_id);
                    genome = address[0];  
                    if (locations[mate][genome] != null && frequenceis[genome] < 20 * Math.log10(genome_sizes[genome])){// should map against this genome
                        sequence = address[1];
                    // calculate the locations based on the offsets in the node    
                        location_array = (int[])r.getProperty(seq_id);
                        seq_len = sequence_length[genome][sequence];
                        if (side == 'F'){
                            for (i = 0; i <= location_array.length - 1; i += 1){
                                if (pointer.canonical ^ is_canonical){
                                    loc = location_array[i] + offset - read_len[mate] + position + K;
                                    if (loc >= 0 && loc <= seq_len - read_len[mate]){
                                        locations[mate][genome].add(new int[]{sequence, -(1 + loc)});
                                    }
                                    //System.out.println("F-" + loc);
                                } else {
                                    loc = location_array[i] + offset - position;
                                    if (loc >= 0 && loc <= seq_len - read_len[mate]){
                                        locations[mate][genome].add(new int[]{sequence, loc});
                                    }
                                    //System.out.println("F+" + loc);
                                }
                            }
                        }else{
                            node_len = (int)node.getProperty("length");
                            for (i = 0; i <= location_array.length - 1; i += 1){
                                if (pointer.canonical ^ is_canonical){
                                    loc = location_array[i] + node_len - K - offset - position;
                                    if (loc >= 0 && loc <= seq_len - read_len[mate]){
                                        locations[mate][genome].add(new int[]{sequence, loc});
                                    }
                                    //System.out.println("R+" + loc);
                                } else {
                                    loc = location_array[i] + node_len - offset - read_len[mate] + position;
                                    if (loc >= 0 && loc <= seq_len - read_len[mate]){
                                        locations[mate][genome].add(new int[]{sequence, -(1 + loc)});
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

        public void cluster_and_examine_locations(int mate, int genome, int sholder){
            int sequence, prev_sequence, prev_start;
            int[] intpair;
            int start, j, k, n, m, count;
            if (locations[mate][genome].size() > 0){
                locations[mate][genome].sort(int_pair_comp);
                n = locations[mate][genome].size();
                for (j = 0; j < n;){
                    intpair = locations[mate][genome].get(j);
                    prev_sequence = sequence = intpair[0];
                    prev_start = intpair[1];
                    for (count = 0; j < n; ++j){
                        intpair = locations[mate][genome].get(j);
                        sequence = intpair[0];
                        start = intpair[1];
                        if (sequence == prev_sequence){
                            if (start - prev_start > sholder){ 
                                hit_counts.add(new int[]{count, prev_start});
                                count = 1;
                                prev_start = start;
                            } else
                                ++count;
                        } else
                            break;
                        //System.out.print(sequence+"_"+start+" ");
                    }
                    //System.out.println("count: " + count);
                    hit_counts.add(new int[]{count, prev_start});
                    hit_counts.sort(int_pair_comp);
                    m = Math.min(hit_counts.size(), MAX_NUM_LOCATIONS);
                    for (k = 0; k < m; ++k){
                        //System.out.print(hit_counts.get(k)[0] + " ");
                        if (sholder == 0)
                            examine_location(mate, genome, prev_sequence, hit_counts.get(k)[1]);
                        else
                            exhaustive_examine_location(mate, genome, prev_sequence, hit_counts.get(k)[1]);
                    }
                    //System.out.println();
                    hit_counts.clear();
                }
                locations[mate][genome].clear();
            }
        }
        
        public void examine_location(int mate, int genome, int sequence, int ref_start){
            boolean forward = ref_start >= 0;
            boolean banded_alignment;
            single_hit h;
            int start, stop;
            if (!forward)
                ref_start = -ref_start - 1;
            start = ref_start - ALIGNMENT_BOUND;
            stop = start + read_len[mate] + 2 * ALIGNMENT_BOUND - 1;
            if (start >= 0 && stop <= sequence_length[genome][sequence] - 1){
                banded_alignment = true;
            } else if(ref_start >= 0 && ref_start + read_len[mate] <= sequence_length[genome][sequence]){
                start = ref_start;
                stop = start + read_len[mate] - 1;
                banded_alignment = false;
            } else
                return;
            num_hits++;
            reference.setLength(0);
            genomeSc.get_sub_sequence(reference, genome, sequence, start, stop - start + 1, true);
            if (alignments[mate].size() < 2 * ALIGNMENT_BOUND * read_len[mate] && 
                find_similar_subject(mate, genome, sequence, start, forward)){
                if (valid_hit())
                    hits[mate][genome].offer(new single_hit(alignment_result));
            } else {
                num_alns++;
                perform_alignment(banded_alignment, mate, genome, sequence, start, forward);
                h = new single_hit(alignment_result);
                alignments[mate].add(h);
                if (valid_hit())
                    hits[mate][genome].offer(h);
            }
        }
        
        public void exhaustive_examine_location(int mate, int genome, int sequence, int ref_start){
            boolean forward = ref_start >= 0;
            single_hit h;
            int start, stop;
            if (!forward)
                ref_start = -ref_start - 1;
            start = Math.max(ref_start - SHOULDER - read_len[mate], 0);
            stop = Math.min(start + 2 * (read_len[mate] + SHOULDER), (int)sequence_length[genome][sequence] - 1);
            num_hits++;
            reference.setLength(0);
            genomeSc.get_sub_sequence(reference, genome, sequence, start, stop - start + 1, true);
            num_alns++;
            perform_alignment(false, mate, genome, sequence, start, forward);
            h = new single_hit(alignment_result);
            if (valid_hit()){
                hits[mate][genome].offer(h);
                /*if (genome == 187){
                    System.out.println(aligner.get_alignment());
                    System.out.println(aligner.get_identity());
                }*/
            }
        }

        public void perform_alignment(boolean banded_alignment, int mate, int genome, int sequence, int start, boolean forward){
            if (banded_alignment){
                bounded_aligner.align(forward?forward_read[mate]:reverse_read[mate], reference); 
                alignment_result.genome = genome;
                alignment_result.sequence = sequence;
                alignment_result.cigar = bounded_aligner.get_cigar().toString();
                alignment_result.identity = bounded_aligner.get_identity();
                alignment_result.score = bounded_aligner.get_similarity();
                alignment_result.start = start;
                alignment_result.offset = bounded_aligner.get_offset();
                alignment_result.length = bounded_aligner.get_range_length();
                alignment_result.deletions = bounded_aligner.get_deletions();
                alignment_result.forward = forward;
                alignment_result.reference = reference.toString();
            } else {
                aligner.align(forward?forward_read[mate]:reverse_read[mate], reference); 
                alignment_result.genome = genome;
                alignment_result.sequence = sequence;
                alignment_result.cigar = aligner.get_cigar().toString();
                alignment_result.identity = aligner.get_identity();
                alignment_result.score = aligner.get_similarity();
                alignment_result.start = start;
                alignment_result.offset = aligner.get_offset();
                alignment_result.length = aligner.get_range_length();
                alignment_result.deletions = aligner.get_deletions();
                alignment_result.forward = forward;
                alignment_result.reference = reference.toString();
            }
        }

        public boolean find_similar_subject(int mate, int curr_genome, int cur_seq, int cur_start, boolean cur_fwd) {
            single_hit similar_alignment;
            boolean found = false;
            Iterator<single_hit> itr = alignments[mate].iterator();
            while (itr.hasNext() && !found) {
                similar_alignment = itr.next();
                if (are_equal(similar_alignment.reference, reference)){
                    alignment_result.genome = curr_genome;
                    alignment_result.sequence = cur_seq;
                    alignment_result.identity = similar_alignment.identity;
                    alignment_result.score = similar_alignment.score;
                    alignment_result.start = cur_start;
                    alignment_result.offset = similar_alignment.offset;
                    alignment_result.length = similar_alignment.length;
                    alignment_result.deletions = similar_alignment.deletions;
                    alignment_result.forward = cur_fwd;
                    alignment_result.cigar = similar_alignment.cigar.toString();
                    alignment_result.reference = reference.toString();
                    found = true;
                }
            }
            return found;
        }  
        
        boolean valid_hit(){
                return (alignment_result.identity > MIN_IDENTITY &&
                        alignment_result.length >= MIN_HIT_LENGTH && 
                        alignment_result.start + alignment_result.offset >= 0 && 
                        alignment_result.start + alignment_result.offset + 
                        alignment_result.deletions + read_len[0] 
                        <= sequence_length[alignment_result.genome][alignment_result.sequence]);            
        }
        
        public boolean are_equal(String s1, StringBuilder s2){
            boolean are_equal = s1.length() == s2.length();
            for (int i = 0; are_equal && i < s1.length(); ++i)
                if (s1.charAt(i) != s2.charAt(i))
                    are_equal = false;
            return are_equal;
        }
        
        public void report_hits(){
            int genome, i;
            if (ALIGNMENT_MODE < 0){ // pan-genomic best
                for (i = 0; i < genome_numbers.size(); ++i){
                    genome = genome_numbers.get(i);
                    collect_hits(genome);
                }
                call_mode();
                clear_hits_list();
            } else {
                //best_hit[0] = null;
                //best_hit[1] = null;
                for (i = 0; i < genome_numbers.size(); ++i){
                    genome = genome_numbers.get(i);
                    collect_hits(genome);
                    call_mode();
                    clear_hits_list();
                }
                //if (unmapped_genomes[0].size() > 0)
                //    check_unmapped_genomes();
            }
            alignments[0].clear();
            if (paired_end)
                alignments[1].clear();
        }
        
        public void check_unmapped_genomes(){
            int i, mate;
            int genome1, genome2;
            single_hit s;
            Relationship rel;
            Node node, neighbor;
            IndexPointer pointer;
            int[] loc2;
            int loc1, seq1, seq2, node_len;
            long high;
            char side1, side2;
            String origin1;
            Iterator<Integer> itr;
            for (mate = 0; mate < num_segments; ++mate){
                s = best_hit[mate];
                if (s != null){
                    genome1 = s.genome;
                    seq1 = s.sequence;
                    try (Transaction tx = graphDb.beginTx()) {
                        loc1 = Math.max(s.start - SHOULDER, 0);
                        pointer = locate(graphDb, genomeSc, indexSc, genome1, seq1, loc1 + 1);
                        origin1 = "G" + genome1 + "S" + seq1;
                        node = graphDb.getNodeById(pointer.node_id);
                        side1 = pointer.canonical ? 'F' : 'R';
                        high = Math.min(s.start + read_len[mate] + SHOULDER - 1, sequence_length[genome1][seq1] - 1);
                        while (loc1 < high){
                            //System.out.println(loc1);
                            node_len = (int) node.getProperty("length");
                            for (Relationship r: node.getRelationships(Direction.INCOMING)){
                                side2 = r.getType().name().charAt(1);
                                for (String origin2: r.getPropertyKeys()){
                                    loc2 = (int[])r.getProperty(origin2);
                                    genome2 = Integer.parseInt(origin2.split("S")[0].substring(1));
                                    //System.out.println(unmapped_genomes[0].size() + " "+genome2);
                                    if (unmapped_genomes[mate].contains(genome2)){
                                        seq2 = Integer.parseInt(origin2.split("S")[1]);
                                        if (side1 == side2){
                                            for (i = 0; i < loc2.length; ++i)
                                                locations[mate][genome2].add(new int[]{seq2, loc2[i]});
                                        } else {
                                            for (i = 0; i < loc2.length; ++i)
                                                locations[mate][genome2].add(new int[]{seq2, -loc2[i]});
                                        }
                                    }
                                }
                            }
                            loc1 += node_len - K + 1;
                            rel = get_outgoing_edge(node, origin1, loc1);
                            if (rel == null)
                                break;
                            else
                                neighbor = rel.getEndNode();
                            node = neighbor;
                            side1 = rel.getType().name().charAt(1);
                        } // while
                        itr = unmapped_genomes[mate].iterator();
                        while (itr.hasNext())
                            cluster_and_examine_locations(mate, itr.next(), SHOULDER + read_len[mate]);
                        tx.success();
                    }
                }
            }
            unmapped_genomes[0].clear();
            if (paired_end)
                unmapped_genomes[1].clear();
        }
        
        public void call_mode(){
            switch (Math.abs(ALIGNMENT_MODE)){
                case 0: // all hits
                    report_all_hit(false);
                break;
                case 1: // unique best hits
                    report_unique_hit();
                break;    
                case 2: // one best hits
                    report_one_hit();
                break;    
                case 3: // all best hits
                    report_all_hit(true);
            }            
        }
        
        public void collect_hits(int genome){
            if (!paired_end){
                if (hits[0][genome].isEmpty()){
                    //unmapped_genomes[0].add(genome);
                    single_hits.offer(new single_hit(genome,0,0,-1,-1,0,0,0,true,null,null));
                } else {
                    while(!hits[0][genome].isEmpty())
                        single_hits.offer(hits[0][genome].remove());
                    /*single_hit h = single_hits.peek();
                    if (best_hit[0] == null || h.score > best_hit[0].score)
                       best_hit[0] = h;*/
                }
            } else {
                boolean reads_paired = false;
                single_hit h1, h2;
                Iterator<single_hit> itr1;
                Iterator<single_hit> itr2;
                int frag_len, best_frag_len;
                if (hits[0][genome].isEmpty() && hits[1][genome].isEmpty()){
                    //unmapped_genomes[0].add(genome);
                    //unmapped_genomes[1].add(genome);
                    paired_hits.offer(new paired_hit(Integer.MAX_VALUE, new single_hit(genome,0,0,-1,-1,0,0,0,true,null,null),new single_hit(genome,0,0,-1,-1,0,0,0,true,null,null)));
                } else if (hits[0][genome].isEmpty() && !hits[1][genome].isEmpty()){
                    //unmapped_genomes[0].add(genome);
                    while (!hits[1][genome].isEmpty())
                        paired_hits.offer(new paired_hit(Integer.MAX_VALUE, new single_hit(genome,0,0,-1,-1,0,0,0,true,null,null),new single_hit(hits[1][genome].remove())));
                    /*h2 = paired_hits.peek().h2;
                    if (best_hit[1] == null || h2.score > best_hit[1].score)
                       best_hit[1] = h2;*/
                } else if (!hits[0][genome].isEmpty() && hits[1][genome].isEmpty()){
                    //unmapped_genomes[1].add(genome);
                    while (!hits[0][genome].isEmpty())
                        paired_hits.offer(new paired_hit(Integer.MAX_VALUE, new single_hit(hits[0][genome].remove()), new single_hit(genome,0,0,-1,-1,0,0,0,true,null,null)));
                    /*h1 = paired_hits.peek().h1;
                    if (best_hit[0] == null || h1.score > best_hit[0].score)
                       best_hit[0] = h1;*/
                } else {
                    itr1 = hits[0][genome].iterator();
                    while (itr1.hasNext()){
                        frag_len = best_frag_len = Integer.MAX_VALUE;
                        h1 = itr1.next();
                        itr2 = hits[1][genome].iterator();
                        while (itr2.hasNext()){
                            h2 = itr2.next();
                            if (h1.sequence != h2.sequence)
                                continue;
                            frag_len = fragment_length(h1, h2);
                            if (frag_len < best_frag_len)
                                best_frag_len = frag_len;
                            else
                                continue;
                            paired_hits.offer(new paired_hit(frag_len, new single_hit(h1), new single_hit(h2)));
                            reads_paired = true;
                        }
                    }
                    if (!reads_paired)
                        paired_hits.offer(new paired_hit(Integer.MAX_VALUE, new single_hit(hits[0][genome].peek()), new single_hit(hits[1][genome].peek())));
                    /*h1 = paired_hits.peek().h1;
                    if (best_hit[0] == null || h1.score > best_hit[0].score)
                       best_hit[0] = h1;
                    h2 = paired_hits.peek().h2;
                    if (best_hit[1] == null || h2.score > best_hit[1].score)
                       best_hit[1] = h2;*/
                    hits[0][genome].clear();
                    hits[1][genome].clear();
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
                        if (shc.compare(h, best_hit) > 0){
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
                best_hit = new paired_hit(h.fragment_length, h.h1, h.h2);
                if (best_hit.get_max_start() != -1){
                    if(!paired_hits.isEmpty()){
                        h = paired_hits.remove();
                        if (phc.compare(h, best_hit) > 0){
                            if (best_hit.get_min_start() != -1){
                                write_paired_sam_record(best_hit, new int[]{1, 1});
                                unique[best_hit.h1.genome]++;
                                unique[best_hit.h2.genome]++;
                            } else if (best_hit.h1.start != -1){
                                write_paired_sam_record(best_hit, new int[]{1, 5});
                                unique[best_hit.h1.genome]++;
                                unmapped[best_hit.h2.genome]++;
                            } else if (best_hit.h2.start != -1){
                                write_paired_sam_record(best_hit, new int[]{5, 1});
                                unique[best_hit.h2.genome]++;
                                unmapped[best_hit.h1.genome]++;
                            }
                        } else {
                            if (best_hit.get_min_start() != -1){
                                write_paired_sam_record(best_hit, new int[]{1, 1});
                                shared[best_hit.h1.genome]++;
                                shared[best_hit.h2.genome]++;
                            } else if (best_hit.h1.start != -1){
                                write_paired_sam_record(best_hit, new int[]{1, 5});
                                shared[best_hit.h1.genome]++;
                                unmapped[best_hit.h2.genome]++;
                            } else if (best_hit.h2.start != -1){
                                write_paired_sam_record(best_hit, new int[]{5, 1});
                                shared[best_hit.h2.genome]++;
                                unmapped[best_hit.h1.genome]++;
                            }
                        }
                    } else {
                        if (best_hit.get_min_start() != -1){
                            write_paired_sam_record(best_hit, new int[]{1, 1});
                            unique[best_hit.h1.genome]++;
                            unique[best_hit.h2.genome]++;
                        } else if (best_hit.h1.start != -1){
                            write_paired_sam_record(best_hit, new int[]{1, 5});
                            unique[best_hit.h1.genome]++;
                            unmapped[best_hit.h2.genome]++;
                        } else if (best_hit.h2.start != -1){
                            write_paired_sam_record(best_hit, new int[]{5, 1});
                            unique[best_hit.h2.genome]++;
                            unmapped[best_hit.h1.genome]++;
                        }
                    }
                } else {
                    unmapped[best_hit.h1.genome]++;
                    unmapped[best_hit.h2.genome]++;
                    write_paired_sam_record(best_hit, new int[]{5, 5});
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
                    if (h.start == -1 || shc.compare(h, best_hit) == 1)
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
                    //if (h.get_identity() == best_hit.get_identity())
                    //System.out.print(h.fragment_length+" ");
                    if (h.get_max_start() == -1 || phc.compare(h, best_hit) == 1)
                        break;
                    sum_freq += raw_abundance[h.h1.genome] / genome_sizes[h.h1.genome];
                    paired_hits_2.add(h);
                }
                //System.out.println();
                rnd = rand.nextDouble();
                while(!paired_hits_2.isEmpty()){
                    best_hit = paired_hits_2.remove();
                    freq += raw_abundance[best_hit.h1.genome] / genome_sizes[best_hit.h1.genome];
                    if (rnd < freq / sum_freq)
                        break;
                }
                if (best_hit.get_max_start() != -1){
                    if (best_hit.get_min_start() != -1){
                        write_paired_sam_record(best_hit, new int[]{1, 1});
                        if (count == 1){
                            unique[best_hit.h1.genome]++;
                            unique[best_hit.h2.genome]++;
                        } else {
                            shared[best_hit.h1.genome]++;
                            shared[best_hit.h2.genome]++;
                        }
                    } else if (best_hit.h1.start != -1){
                        write_paired_sam_record(best_hit, new int[]{1, 5});
                        unmapped[best_hit.h2.genome]++;
                        if (count == 1)
                            unique[best_hit.h1.genome]++;
                        else 
                            shared[best_hit.h1.genome]++;
                    } else if (best_hit.h2.start != -1){
                        write_paired_sam_record(best_hit, new int[]{5, 1});
                        unmapped[best_hit.h1.genome]++;
                        if (count == 1)
                            unique[best_hit.h2.genome]++;
                        else 
                            shared[best_hit.h2.genome]++;
                    }
                } else { 
                       write_paired_sam_record(best_hit, new int[]{5, 5});
                        unmapped[best_hit.h1.genome]++;
                        unmapped[best_hit.h2.genome]++;
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
                    if (h.start == -1 || (best && shc.compare(h, best_hit) > 0))
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
                    if (h.get_max_start() == -1 || (best && phc.compare(h, best_hit) > 0))
                        break;
                    paired_hits_2.add(h);
                }
                if (best_hit.get_max_start() != -1){
                    prev_genome = -1;
                    while (!paired_hits_2.isEmpty()){
                        best_hit = paired_hits_2.remove();
                        if (best_hit.get_min_start() != -1){
                            write_paired_sam_record(best_hit, new int[]{best_hit.h1.genome == prev_genome ? 257 : 1, best_hit.h2.genome == prev_genome ? 257 : 1});
                            if (count == 1){
                                unique[best_hit.h1.genome]++;
                                unique[best_hit.h2.genome]++;
                            } else {
                                shared[best_hit.h1.genome]++;
                                shared[best_hit.h2.genome]++;
                            }
                        } else if (best_hit.h1.start != -1){
                            write_paired_sam_record(best_hit, new int[]{best_hit.h1.genome == prev_genome ? 257 : 1, 5});
                            unmapped[best_hit.h2.genome]++;
                            if (count == 1)
                                unique[best_hit.h1.genome]++;
                            else 
                                shared[best_hit.h1.genome]++;
                        } else if (best_hit.h2.start != -1){
                            write_paired_sam_record(best_hit, new int[]{5, best_hit.h2.genome == prev_genome ? 257 : 1});
                            unmapped[best_hit.h1.genome]++;
                            if (count == 1)
                                unique[best_hit.h2.genome]++;
                            else 
                                shared[best_hit.h2.genome]++;
                        }
                        prev_genome = best_hit.h1.genome;
                    }
                } else { 
                       write_paired_sam_record(best_hit, new int[]{5, 5});
                        unmapped[best_hit.h1.genome]++;
                        unmapped[best_hit.h2.genome]++;
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

        public void write_single_sam_record(single_hit h, int flag){
            SAMRecord sam_record = new SAMRecord(null);
            if (h.start == -1){
                sam_record.setReadName(reads[0].getReadName());
                sam_record.setFlags(flag);
                sam_record.setReadString(forward_read[0].toString());
            } else {
                flag |= h.forward?0:16; // direction
                sam_record.setReferenceName(sequence_titles[h.genome][h.sequence].split("\\s")[0]);
                sam_record.setFlags(flag);
                sam_record.setReadName(reads[0].getReadName());
                sam_record.setAlignmentStart(h.start + h.offset + 1);
                sam_record.setMappingQuality((int)Math.ceil(-10 * Math.log10(1 - (h.identity - 0.001))));
                sam_record.setCigarString(h.cigar);
                sam_record.setReadString((h.forward?forward_read[0]:reverse_read[0]).toString());
                sam_record.setBaseQualityString(quality[0].toString());
            }
            //synchronized(sam_headers[h.genome]){
                synchronized(sam_writers[h.genome]){
                    //sam_record.setHeader(sam_headers[h.genome]);
                    sam_writers[h.genome].addAlignment(sam_record);
                }
            //}
        }
        
        public void write_paired_sam_record(paired_hit h, int[] flag){
        // the first segment    
            SAMRecord sam_record1, sam_record2;
            sam_record1 = new SAMRecord(null);
            int position1 = h.h1.start + h.h1.offset + 1;
            int position2 = h.h2.start + h.h2.offset + 1;;
            flag[0] |= 64; 
            if (h.h1.start == -1){
                sam_record1.setReadName(read_name[0]);
                sam_record1.setFlags(flag[0]);
                sam_record1.setReadString(forward_read[0].toString());
            } else {
                flag[0] |= h.h1.forward?0:16; // SEQ being reverse complemented
                flag[0] |= !h.h2.forward?32:0; // SEQ being reverse complemented
            //qname
                sam_record1.setReadName(read_name[0]); 
            //flag    
                sam_record1.setFlags(flag[0]); 
            //rname    
                sam_record1.setReferenceName(sequence_titles[h.h1.genome][h.h1.sequence].split("\\s")[0]);
            //pos    
                sam_record1.setAlignmentStart(position1);
            //mapq    
                sam_record1.setMappingQuality((int)Math.ceil(-10 * Math.log10(1 - (h.h1.identity - 0.001))));
                //sam_record.setMappingQuality((int)Math.round(h.identity1 * 100));
            //cigar    
                sam_record1.setCigarString(h.h1.cigar);
            //rnext
                if (h.h2.start == -1)
                    sam_record1.setMateReferenceName("*");
                else 
                    sam_record1.setMateReferenceName(sequence_titles[h.h2.genome][h.h2.sequence].split("\\s")[0]);
            //pnext 
                sam_record1.setMateAlignmentStart(position2);
            //tlen
                if (h.h2.start == -1)
                    sam_record1.setInferredInsertSize(0);
                else if (position2 > position1)
                    sam_record1.setInferredInsertSize(position2 - position1 + reads[0].getReadLength());
                else
                    sam_record1.setInferredInsertSize(-(position1 - position2 + reads[1].getReadLength()));
            //seq    
                sam_record1.setReadString((h.h1.forward?forward_read[0]:reverse_read[0]).toString());
            //qual    
                sam_record1.setBaseQualityString(quality[0].toString());
            }
        // the second segment    
            sam_record2 = new SAMRecord(null);
            flag[1] |= 128; 
            if (h.h2.start == -1){ // not mapped
                sam_record2.setReadName(read_name[1]);
                sam_record2.setFlags(flag[1]);
                sam_record2.setReadString(forward_read[1].toString());
            } else {
                flag[1] |= h.h2.forward?0:16; // SEQ being reverse complemented
                flag[1] |= !h.h1.forward?32:0; // SEQ being reverse complemented
            //qname
                sam_record2.setReadName(read_name[1]); 
            //flag    
                sam_record2.setFlags(flag[1]); 
            //rname    
                sam_record2.setReferenceName(sequence_titles[h.h2.genome][h.h2.sequence].split("\\s")[0]);
            //pos    
                sam_record2.setAlignmentStart(position2);
            //mapq    
                sam_record2.setMappingQuality((int)Math.ceil(-10 * Math.log10(1 - (h.h2.identity - 0.001))));
                //sam_record.setMappingQuality((int)Math.round(h.identity2 * 100));
            //cigar    
                sam_record2.setCigarString(h.h2.cigar);
            //rnext
                if (h.h1.start == -1)
                    sam_record2.setMateReferenceName("*");
                else 
                    sam_record2.setMateReferenceName(sequence_titles[h.h1.genome][h.h1.sequence].split("\\s")[0]);
            //pnext    
                sam_record2.setMateAlignmentStart(position1);
            //tlen
                if (h.h1.start == -1)
                    sam_record2.setInferredInsertSize(0);
                else if (position2 > position1)
                    sam_record2.setInferredInsertSize(-(position2 - position1 + reads[0].getReadLength()));
                else
                    sam_record2.setInferredInsertSize(position1 - position2 + reads[1].getReadLength());
            //seq    
                sam_record2.setReadString((h.h2.forward?forward_read[1]:reverse_read[1]).toString());
            //qual    
                sam_record2.setBaseQualityString(quality[1].toString());
            }
            synchronized(sam_writers[h.h1.genome]){
                //sam_record.setHeader(sam_headers[h.h1.genome]);
                sam_writers[h.h1.genome].addAlignment(sam_record1);
                sam_writers[h.h2.genome].addAlignment(sam_record2);
            }
        }
        

        public int fragment_length(single_hit h1, single_hit h2){
            int frag_len;
            int position1, position2;
            position1 = h1.start + h1.offset;
            position2 = h2.start + h2.offset;
            if (h1.forward)
                frag_len = position2 + reads[1].getReadLength() - position1;
            else
                frag_len = position1 + reads[0].getReadLength() - position2;
            if (frag_len < Math.max(read_len[0], read_len[1]) || frag_len > MAX_FRAGMENT_LENGTH )
                frag_len = Integer.MAX_VALUE;
            return frag_len;
        }
    }   
    
    public void map_reads() {
        int i, j, t, number, genome, n = 0;
        Node pangenome_node;
        ArrayList<Integer>[] genome_numbers;
        String line;
        BufferedReader in;
        BufferedWriter out;
        int total_unique, total_mapped, total_unmapped;
        boolean paired;
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
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();//.setConfig(GraphDatabaseSettings.pagecache_memory, "64g")
        registerShutdownHook(graphDb);
        indexDb = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH, "sorted");
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
        number_of_alignments = new AtomicInteger(0);
        number_of_hits = new AtomicInteger(0);
        genome_numbers = new ArrayList[THREADS];
        for (t = 0; t < THREADS; ++t)
            genome_numbers[t] = new ArrayList();
        try {
            in = new BufferedReader(new FileReader(PATH_TO_THE_GENOME_NUMBERS_FILE));
            for (n = 0; in.ready(); ){
                line = in.readLine().trim();
                if (!line.equals("")){
                    number = Integer.parseInt(line);
                    if (number > 0 && number <= genomeDb.num_genomes){
                       for (t = 0; t < THREADS; ++t)
                            genome_numbers[t].add(number);
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
        for (i = 0; i < genome_numbers[0].size(); ++i){
            genome = genome_numbers[0].get(i);
            headers[genome] = new SAMFileHeader();
            for (j = 1; j <= genomeDb.num_sequences[genome]; ++j)
                headers[genome].addSequence(new SAMSequenceRecord(genomeDb.sequence_titles[genome][j].split("\\s")[0], (int)genomeDb.sequence_length[genome][j]));
            headers[genome].addProgramRecord(new SAMProgramRecord("PanTools"));
            if (BAMFORMAT)
                sams[genome] = new SAMFileWriterFactory().makeBAMWriter(headers[genome], false, 
                        new File(OUTPUT_PATH + "/pantools_" + genome + ".bam"));
            else
                sams[genome] = new SAMFileWriterFactory().makeSAMWriter(headers[genome], false, 
                        new File(OUTPUT_PATH + "/pantools_" + genome + ".sam"));
        }
        fastq_records = new LinkedBlockingQueue[THREADS][];
        for (t = 0; t < THREADS; ++t){
            fastq_records[t] = new LinkedBlockingQueue[2];
            fastq_records[t][0] = new LinkedBlockingQueue((int)(heapSize/(20000*THREADS)));
            if (paired)
                fastq_records[t][1] = new LinkedBlockingQueue((int)(heapSize/(20000*THREADS)));
        }
        if (n > 0){
            System.out.println("\nMapping" + (paired?" paired-end ": " ") + "reads on " + n + " genome(s) :");
            try{
                ExecutorService es = Executors.newFixedThreadPool(THREADS + 1);
                es.execute(new Generate_fatsq_records(paired));
                for (t = 0; t < THREADS; ++t)
                    es.execute(new Map(t, genome_numbers[t], paired, sams));
                es.shutdown();
                es.awaitTermination(10, TimeUnit.DAYS);        
            } catch (InterruptedException e){

            }
            

            total_unique = total_mapped = total_unmapped = 0;
            try{
                out = new BufferedWriter(new FileWriter(OUTPUT_PATH + "/mapping_summary.txt"));
                out.write("Genome\tTotal\tUnique\tUnmapped\n");
                System.out.print("\nGenome\tTotal\tUnique\tUnmapped\n");
                for (i = 0; i < genome_numbers[0].size(); ++i){
                    genome = genome_numbers[0].get(i);
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
            System.out.println("........................................");
            System.out.println("\t" + total_mapped + "\t" + total_unique + "\t" + total_unmapped);
            System.out.println("Number_of_hits = " + number_of_hits);
            System.out.println("Number_of_alignments = " + number_of_alignments + "\n");
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
        add_sequence_properties();
        localize_nodes();
        System.out.println("Number of kmers:   " + indexSc.length());
        System.out.println("Number of nodes:   " + num_nodes);
        System.out.println("Number of edges:   " + num_edges);
        System.out.println("Number of bases:   " + num_bases);
        System.out.println("Number of degenerate nodes:   " + num_degenerates);
        try (Transaction tx = graphDb.beginTx()) {
            pangenome_node.setProperty("k_mer_size", K_SIZE);
            pangenome_node.setProperty("num_k_mers", indexSc.length());
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
        indexSc = new IndexScanner(indexDb);
        genomeSc = new SequenceScanner(genomeDb, previous_num_genomes + 1, 1, K_SIZE, indexSc.get_pre_len());
    // the sequences should be dropped out as they will change and add_sequence_properties() function will rebuild them.    
        drop_nodes_property("sequence");
        drop_nodes_property("frequencies");
        drop_nodes_property("frequency");
    // the edge colors should be dropped out as they will change and localize_nodes() function will rebuild them again.    
        drop_edges_colors();
        construct_pangenome(pangenome_node);
        add_sequence_properties();
        localize_nodes();
        System.out.println("Number of kmers:   " + indexSc.length());
        System.out.println("Number of nodes:   " + num_nodes);
        System.out.println("Number of edges:   " + num_edges);
        System.out.println("Number of bases:   " + num_bases);
        System.out.println("Number of degenerate nodes:   " + num_degenerates);
        try (Transaction tx = graphDb.beginTx()) {
            pangenome_node.setProperty("k_mer_size", K_SIZE);
            pangenome_node.setProperty("num_k_mers", indexSc.length());
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
        int genome, seqience, begin, end, j, len;
        long byte_number = 0;
        SequenceDatabase genomeDb = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH, graphDb);
        StringBuilder seq = new StringBuilder();
        for (genome = 1; genome <= genomeDb.num_genomes; ++genome) {
            for (seqience = 1; seqience <= genomeDb.num_sequences[genome]; ++seqience) {
                begin = 1;
                end = (int) genomeDb.sequence_length[genome][seqience];
                extract_sequence_from_graph(seq, genome, seqience, begin, end);
                len = seq.length();
                if (len % 2 == 1) {
                    --len;
                }
                for (j = 0; j < len; j += 2, ++byte_number) {
                    genomeDb.genomes_buff[(int) (byte_number / genomeDb.MAX_BYTE_COUNT)].put((byte) ((binary[seq.charAt(j)] << 4) | binary[seq.charAt(j + 1)]));
                }
                if (len == seq.length() - 1) {
                    genomeDb.genomes_buff[(int) (byte_number / genomeDb.MAX_BYTE_COUNT)].put((byte) (binary[seq.charAt(len)] << 4));
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
                } else {
                    try {
                        FileUtils.deleteRecursively(new File(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH));
                    } catch (IOException ioe) {
                        System.out.println("Failed to delete the index database!");
                        System.exit(1);  
                    }                    
                    indexDb = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH, PATH_TO_THE_GENOMES_FILE, genomeDb, K_SIZE);
                }
            } else
                indexDb = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH, PATH_TO_THE_GENOMES_FILE, genomeDb, K_SIZE);
            indexSc = new IndexScanner(indexDb);
            K_SIZE = indexSc.get_K();
            genomeSc = new SequenceScanner(genomeDb, 1, 1, K_SIZE, indexSc.get_pre_len());
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
                indexSc = new IndexScanner(indexDb);
                K_SIZE = indexSc.get_K();
                genomeSc = new SequenceScanner(genomeDb, 1, 1, K_SIZE, indexSc.get_pre_len());
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
        StringBuilder seq;
        int genome, sequence, begin, end, num_regions = 0, proper_regions = 0;
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
        registerShutdownHook(graphDb);
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
        indexSc = new IndexScanner(indexDb);
        K_SIZE = indexSc.get_K();
        genomeSc = new SequenceScanner(genomeDb, 1, 1, K_SIZE, indexSc.get_pre_len());
        try (Transaction tx = graphDb.beginTx()) {
            try (BufferedReader in = new BufferedReader(new FileReader(PATH_TO_THE_REGIONS_FILE))) {
                fields = PATH_TO_THE_REGIONS_FILE.split("\\/");
                out_file_name = PATH_TO_THE_PANGENOME_DATABASE + "/" + fields[fields.length - 1] + ".fasta";
                BufferedWriter out = new BufferedWriter(new FileWriter(out_file_name));
                while (in.ready()) {
                    line = in.readLine().trim();
                    if (line.equals("")) {
                        continue;
                    }
                    fields = line.trim().split("\\s");
                    genome = Integer.parseInt(fields[0]);
                    sequence = Integer.parseInt(fields[1]);
                    begin = Integer.parseInt(fields[2]);
                    end = Integer.parseInt(fields[3]);
                    if (genome <= genomeDb.num_genomes && sequence <= genomeDb.num_sequences[genome] && begin >= 1 && end <= genomeDb.sequence_length[genome][sequence]){
                        proper_regions++;
                        out.write(">genome:" + genome + " sequence:" + sequence + " from:" + begin + " to:" + end + " length:" + (end - begin + 1) + "\n");
                        begin -= 1;
                        end -= 1;
                        seq.setLength(0);
                        genomeSc.get_sub_sequence(seq, genome, sequence, begin, end, true);
                        write_fasta(out, seq.toString(), 70);
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
        String genome_number;
        int genome, sequence, begin, end;
        StringBuilder seq;
        Node pangenome_node;
        graphDb = new GraphDatabaseFactory().newEmbeddedDatabaseBuilder(new File(PATH_TO_THE_PANGENOME_DATABASE + GRAPH_DATABASE_PATH))
                .setConfig(keep_logical_logs, "4 files").newGraphDatabase();
        registerShutdownHook(graphDb);
        startTime = System.currentTimeMillis();
        genomeDb = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH);
        indexDb = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH, "sorted");
        indexSc = new IndexScanner(indexDb);
        K_SIZE = indexSc.get_K();
        genomeSc = new SequenceScanner(genomeDb, 1, 1, K_SIZE, indexSc.get_pre_len());
        seq = new StringBuilder();
        try (Transaction tx = graphDb.beginTx()) {
            pangenome_node = graphDb.findNodes(pangenome_label).next();
            if (pangenome_node == null) {
                System.out.println("Can not locate database node!");
                System.exit(1);
            }
            try {
                in = new BufferedReader(new FileReader(PATH_TO_THE_GENOME_NUMBERS_FILE));
                while (in.ready()) {
                    genome_number = in.readLine().trim();
                    if (genome_number.equals(""))
                        continue;
                    try{
                        genome = Integer.parseInt(genome_number);
                    }catch(NumberFormatException e){
                        System.out.println(genome_number + "is not a valid genome number.");
                        continue;
                    }
                    if (genome < 1 || genome > genomeDb.num_genomes){
                        System.out.println(genome_number + "is not a valid genome number.");
                        continue;
                    }
                    System.out.println("Reconstructing genome " + genome_number + "...");
                    try {
                        out = new BufferedWriter(new FileWriter(PATH_TO_THE_PANGENOME_DATABASE + "/Genome_" + genome_number + ".fasta"));
                        for (sequence = 1; sequence <= genomeDb.num_sequences[genome]; ++sequence) {
                            System.out.println("Sequence " + sequence + " length = " + genomeDb.sequence_length[genome][sequence]);
                            //begin = 1;
                            //end = (int)genomeDb.sequence_length[genome][sequence];
                            out.write(">" + genomeDb.sequence_titles[genome][sequence] + "\n");
                            begin = 0;
                            end = (int)genomeDb.sequence_length[genome][sequence] - 1;
                            seq.setLength(0);
                            genomeSc.get_sub_sequence(seq, genome, sequence, begin, end, true);
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
        indexDb = new IndexDatabase(PATH_TO_THE_PANGENOME_DATABASE + INDEX_DATABASE_PATH, "sorted");
        indexSc = new IndexScanner(indexDb);
        genomeDb = new SequenceDatabase(PATH_TO_THE_PANGENOME_DATABASE + GENOME_DATABASE_PATH);
        startTime = System.currentTimeMillis();
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
                    start = locate(graphDb, genomeSc, indexSc, genome1, seq1, loc1 + 1);
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
    public void extract_sequence_from_graph(StringBuilder seq, int genome, int sequence, int begin, int end) {
        Relationship rel;
        Node neighbor, node;
        IndexPointer start_ptr;
        int loc, len = 0, node_len, neighbor_len, seq_len, position;
        String rel_name, origin;
        --begin;
        --end;
        origin = "G" + genome + "S" + sequence;
        seq_len = end - begin + 1;
        seq.setLength(0);
        start_ptr = locate(graphDb, genomeSc, indexSc, genome, sequence, begin);
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
            loc = (begin + len) - K_SIZE + 1;
            rel = get_outgoing_edge(node, origin, loc);
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
    public static IndexPointer locate(GraphDatabaseService graphDb, SequenceScanner genomeSc, IndexScanner indexSc, int genome, int sequence, int position) {
        int i, code, node_start_pos, low, high, mid , node_len, genomic_pos, k_size, pre_len;
        boolean forward, degenerate;
        Node node, neighbor, seq_node;
        Relationship rel;
        String anchor_sides, origin;
        long[] anchor_nodes;
        int[] anchor_positions;
        IndexPointer pointer;
        k_size = indexSc.get_K();
        pre_len = indexSc.get_pre_len();
        kmer first_kmer = new kmer(k_size, pre_len);
        degenerate = false;
        for (i = 0; i < k_size && !degenerate; ++i){
            code = genomeSc.get_code(genome, sequence, position + i);
            if (code < 4)
                first_kmer.next_kmer(code);
            else
               degenerate = true; 
        }
        if (!degenerate){
            pointer = new IndexPointer();
            indexSc.get_pointer(pointer, indexSc.find(first_kmer));
            pointer.canonical ^= first_kmer.get_canonical();
            return pointer;
        }
        origin = "G" + genome + "S" + sequence;
        genomic_pos = genome - 1;
        seq_node = graphDb.findNode(sequence_label, "identifier", genome+"_"+sequence);
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
            while (node_start_pos + node_len <= genomic_pos) {
                genome = node_start_pos + node_len - k_size + 1;
                rel = get_outgoing_edge(node, origin, genome);
                if (rel == null){
                    System.out.println("Failed to locate address : " + genome + " " + sequence + " "+ genome);
                    break;
                }
                neighbor = rel.getEndNode();
                forward = rel.getType().name().charAt(1) == 'F';
                node_start_pos += node_len - k_size + 1;
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
        node_last_kmer=indexSc.find(genomeSc.make_kmer(gen,seq,loc+pos-1));
        split_first_kmer=indexSc.get_next_index(node_last_kmer);
        indexSc.put_next_index(-1L, node_last_kmer); 
        split_node.setProperty("first_kmer",split_first_kmer);
        split_node.setProperty("last_kmer",node.getProperty("last_kmer"));
        s_id=(int)split_node.getId();
        for(i=0,inx=split_first_kmer;inx!=-1L;inx=indexSc.get_next_index(inx),++i) // update kmer coordinates
        {
            indexSc.put_node_id(s_id, inx);
            indexSc.put_position(i, inx);
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
        address = genomeSc.get_address();
        address[2] -= K_SIZE - 1;
        ++num_nodes;
        Node new_node = graphDb.createNode(nucleotide_label);
        node_id = new_node.getId();
        if (DEBUG) System.out.println("create "+new_node.getId());
        new_node.setProperty("address", address);
        new_node.setProperty("length", K_SIZE);
        new_node.setProperty("last_kmer",genomeSc.get_curr_index());
        new_node.setProperty("first_kmer",genomeSc.get_curr_index());
    // Set the pointer to the Kmer in the pointer database    
        indexSc.put_pointer(new_node.getId(), 0, genomeSc.get_curr_kmer().get_canonical(), -1l, genomeSc.get_curr_index());
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
        while (genomeSc.get_position() < genomeSc.get_sequence_length() - 1) { // Not reached to the end of the sequence
            if (DEBUG) System.out.println("extend " + genomeSc.get_position());
            if (genomeSc.get_code(1) > 3) { // hit a degenerate region
                genomeSc.next_position();
                begin = genomeSc.get_position() - K_SIZE + 1;
                curr_node.setProperty("length", len);
                curr_node.setProperty("last_kmer",last_kmer);                
                genomeSc.jump_forward();
                if (genomeSc.get_position() >= genomeSc.get_sequence_length() - 1){
                    genomeSc.next_position();// to acheive the right length for the degenerate node    
                    finish = true;
                }                
                int[] add = genomeSc.get_address();
                add[2] = begin;
                degenerate_node = create_degenerate(add);
                connect(curr_node ,degenerate_node, RelTypes.FF); 
                ++num_edges;
                curr_node = degenerate_node;
                break;
            } else {
                genomeSc.next_position();
                genomeSc.get_curr_kmer().next_kmer(genomeSc.get_code(0));
                if (SHOW_KMERS) System.out.println(genomeSc.get_curr_kmer().toString());
            }
            genomeSc.set_curr_index(indexSc.find(genomeSc.get_curr_kmer()));
            if (indexSc.get_node_id(genomeSc.get_curr_index()) == -1L) {
                indexSc.put_next_index(genomeSc.get_curr_index(),last_kmer);
                ++len;
                indexSc.put_pointer(node_id, len - K_SIZE, genomeSc.get_curr_kmer().get_canonical(), -1l, genomeSc.get_curr_index());
                last_kmer=genomeSc.get_curr_index();
            } else {
                broke = true;
                break;
            }
        }
        if (degenerate_node == null){
            new_node.setProperty("length", len);
            new_node.setProperty("last_kmer",last_kmer);
        }
        if (!broke && genomeSc.get_position() >= genomeSc.get_sequence_length() - 1) {// because the last kmer is somewhere in the graph we should get connected to
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
        for (pos = 0; pos <= l && genomeSc.get_position() <= genomeSc.get_sequence_length() - 1 && genomeSc.get_code(g, s, loc + pos + K_SIZE - 1) == genomeSc.get_code(0); ++pos) {
            genomeSc.next_position();
         // If hit a degenarate region aplit and branch to a degenerate node 
            if (genomeSc.get_position() <= genomeSc.get_sequence_length() - 1 && genomeSc.get_code(0) > 3) {
                begin = genomeSc.get_position() - K_SIZE + 1;
                genomeSc.jump_forward();
                if (genomeSc.get_position() >= genomeSc.get_sequence_length() - 1){
                    genomeSc.next_position();// to acheive the right length for the degenerate node    
                    finish = true;
                }                
                if (pos + 1 <= l) {
                    split_node2 = split(node, pos + 1);
                    if (loop)
                        src = split_node2;
                }
                int[] add = genomeSc.get_address();
                add[2] = begin;
                degenerate_node = create_degenerate(add);
                connect(node ,degenerate_node, RelTypes.FF);
                ++num_edges;
                break;
            }
        }
        if (genomeSc.get_position() == genomeSc.get_sequence_length()) {
            finish = true;
        } else if (degenerate_node == null) // build the Kmer of difference 
            initialize(genomeSc.get_position() - K_SIZE + 1);
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
        for (pos = (int) node.getProperty("length") - K_SIZE; pos >= 0 && genomeSc.get_position() <= genomeSc.get_sequence_length() - 1 && genomeSc.get_code(g, s, loc + pos) == genomeSc.get_complement_current_code(0); --pos) {
            genomeSc.next_position();
            if (genomeSc.get_position() <= genomeSc.get_sequence_length() - 1 && genomeSc.get_code(0) > 3) {
                begin = genomeSc.get_position() - K_SIZE + 1;
                genomeSc.jump_forward();
                if (genomeSc.get_position() >= genomeSc.get_sequence_length() - 1){
                    genomeSc.next_position();// to acheive the right length for the degenerate node    
                    finish = true;
                }                
                if (pos > 0) {
                    split_node2 = split(node, pos);
                    des = split_node2;
                    if (!first_split && loop)
                        src = split_node2;
                }
                int[] add = genomeSc.get_address();
                add[2] = begin;
                degenerate_node = create_degenerate(add);
                connect(split_node2, degenerate_node, RelTypes.RF);
                ++num_edges;
                break;
            }
        }
        if (genomeSc.get_position() == genomeSc.get_sequence_length()) {
            finish = true;
        } else if (degenerate_node == null) // build the Kmer of difference 
            initialize(genomeSc.get_position() - K_SIZE + 1);
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
        if (DEBUG) System.out.println("create_degenerate:"+degenerate_node.getId()+" position:"+genomeSc.get_position()+" begin:"+address[2]);
        degenerate_node.setProperty("address", address);
        degenerate_node.setProperty("length", genomeSc.get_position() - address[2]);
        return degenerate_node;
    }
    
    private void initialize(int start){
        Node degenerate_node;
        if (genomeSc.initialize_kmer(start) < 0){// start with a degenerate kmer
            if (genomeSc.get_position() >= genomeSc.get_sequence_length() - 1){
                genomeSc.next_position();// to acheive the right length for the degenerate node    
                finish = true;
            }                
            int[] add = genomeSc.get_address();
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
            while (!genomeSc.end_of_scan()) {
                System.out.println("Processing genome " + genomeSc.get_genome() + " (" + genomeDb.genome_names[genomeSc.get_genome()] +") :");
                genome_node = graphDb.createNode(genome_label);
                genome_node.setProperty("path", genomeDb.genome_names[genomeSc.get_genome()]);
                genome_node.setProperty("number", genomeSc.get_genome());
                genome_node.setProperty("num_sequences", genomeDb.num_sequences[genomeSc.get_genome()]);
                genome_node.setProperty("date", new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(new Date()));
                pangenome_node.createRelationshipTo(genome_node, RelTypes.has);
                while (!genomeSc.end_of_genome()) {
                    sequence_node = graphDb.createNode(sequence_label);
                    sequence_node.setProperty("genome", genomeSc.get_genome());
                    sequence_node.setProperty("number", genomeSc.get_sequence());
                    sequence_node.setProperty("identifier", genomeSc.get_genome() + "_" + genomeSc.get_sequence());
                    sequence_node.setProperty("title", genomeSc.get_title());
                    sequence_node.setProperty("length", genomeSc.get_sequence_length());
                    sequence_node.setProperty("offset", genomeSc.get_offset());
                    genome_node.createRelationshipTo(sequence_node, RelTypes.has);
                    finish = false;
                    System.out.println("sequence " + genomeSc.get_sequence() + "/" + genomeDb.num_sequences[genomeSc.get_genome()] + 
                    " of genome " + genomeSc.get_genome() + "\tlength=" + genomeSc.get_sequence_length());
                    curr_node = sequence_node;
                    curr_side = 0;
                    if (genomeSc.get_sequence_length() >= K_SIZE){
                        initialize(0);
                        while (!finish) {
                            genomeSc.set_curr_index(indexSc.find(genomeSc.get_curr_kmer()));
                            indexSc.get_pointer(pointer, genomeSc.get_curr_index());
                            if (pointer.node_id == -1L) // kmer is new
                                create_extend();
                            else if (genomeSc.get_curr_kmer().get_canonical() ^ pointer.canonical)// if sides don't agree
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
                    genomeSc.next_sequence();
                }//sequences
                System.out.println((System.currentTimeMillis() - phaseTime) / 1000 + " seconds elapsed.");
                genomeSc.next_genome();
            }//genomes
            tx.success();
        } finally {
            tx.close();
        }
    }
    
    /**
     * To add list of anchor nodes, anchor sides and anchor positions to each sequence_node.
     * These are used for locating genomic regions very quickly. 
     */
    void localize_nodes() {
        ResourceIterator<Node> sequence_iterator;
        LinkedList<Node> sequence_nodes;
        int trsc = 0, i, len, m, neighbor_length = 0;
        long[] frequencies;
        long frequency;
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
                            if ((is_node && genomeSc.compare(address, addr, K_SIZE - 1,
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
                                } else {
                                    initial_coordinate[0] = address[2];
                                    r.setProperty(origin, initial_coordinate);
                                }
                                frequencies = (long[])neighbor.getProperty("frequencies", new long[genomeDb.num_genomes + 1]);
                                ++frequencies[address[0]];
                                neighbor.setProperty("frequencies", frequencies);
                                frequency = (long)neighbor.getProperty("frequency", 0l);
                                ++frequency;
                                neighbor.setProperty("frequency", frequency);
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
                genomeSc.get_sub_sequence(sequence, addr[0], addr[1], addr[2], node_length, true);
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
