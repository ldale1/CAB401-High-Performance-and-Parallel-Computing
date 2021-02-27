package promoters;

import edu.au.jacobi.pattern.Match;
import edu.au.jacobi.pattern.Series;
import qut.*;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.locks.ReentrantLock;

public class Parallel_v2 extends Parallel {

    // Threading
    protected HashMap<String, Sigma70Consensus> consensus = new HashMap<>();
    public ThreadLocal<Series> sigma70_pattern = ThreadLocal.withInitial(() -> Sigma70Definition.getSeriesAll_Unanchored(0.7));
    public final ReentrantLock lock = new ReentrantLock();

    // Constructor
    public Parallel_v2(String referenceFile, String dir) {
        super(referenceFile, dir);
    }

    @Override
    protected List<Gene> ParseReferenceGenes(String referenceFile) throws IOException {
        // Read from an input file
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(referenceFile)));
        List<Gene> referenceGenes = new ArrayList<Gene>();
        while (true) {
            String name = reader.readLine();
            if (name == null)
                break;
            String sequence = reader.readLine();
            referenceGenes.add(new Gene(name, 0, 0, sequence));
            consensus.put(name, new Sigma70Consensus());
        }
        consensus.put("all", new Sigma70Consensus());
        reader.close();
        return referenceGenes;
    }

    @Override
    public void DisplayResults() {
        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet())
            System.out.println(entry.getKey() + " " + entry.getValue());
    }

    @Override
    public String ResultsString() {
        return consensus.values().toString();
    }

    @Override
    public Match PredictPromoter(NucleotideSequence upStreamRegion) {
        return BioPatterns.getBestMatch(sigma70_pattern.get(), upStreamRegion.toString());
    }

    private class GeneThreadFine implements Callable<Void> {
        private Gene gene;
        private Gene referenceGene;
        private NucleotideSequence nucleotides;

        // Constructor
        public GeneThreadFine(Gene gene, Gene referenceGene, NucleotideSequence nucleotides) {
            this.gene = gene;
            this.referenceGene = referenceGene;
            this.nucleotides = nucleotides;
        }

        // Invoke
        @Override
        public Void call() {
            // Copy paste of inner for each loop, with locking for shared heap object
            if (Homologous(gene.sequence, referenceGene.sequence)) {
                Match prediction = PredictPromoter(GetUpstreamRegion(nucleotides, gene));
                if (prediction != null) {
                    lock.lock();
                    try {
                        consensus.get(referenceGene.name).addMatch(prediction);
                        consensus.get("all").addMatch(prediction);
                    }
                    finally {
                        lock.unlock();
                    }
                }
            }
            return null;
        }
    }

    @Override
    public void Run(Integer threads) throws IOException {
        // Scheduling setup
        ExecutorService executorService = Executors.newFixedThreadPool(threads);
        List<Callable<Void>> callableList = new ArrayList<>();

        // Genes in reference file
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);

        // For each Ecoli file
        for (String filename : ListGenbankFiles(dir)) {
            GenbankRecord record = Parse(filename);
            // For each gene in the reference file
            for (Gene referenceGene : referenceGenes) {
                // For each gene in the Ecoli file
                for (Gene gene : record.genes) {
                    callableList.add(new GeneThreadFine(gene, referenceGene, record.nucleotides));
                }
            }
        }

        // Run all scheduled tasks
        try {
            executorService.invokeAll(callableList);
        } catch (InterruptedException ex) {
            ex.printStackTrace();
        }
        executorService.shutdown();
    }
}