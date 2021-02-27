package promoters;

import edu.au.jacobi.pattern.Match;
import edu.au.jacobi.pattern.Series;
import qut.*;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class Parallel_v3 extends Parallel {

    // Threading
    public ConcurrentHashMap<String, Sigma70Consensus> consensus = new ConcurrentHashMap<>();
    public ThreadLocal<Series> sigma70_pattern = ThreadLocal.withInitial(() -> Sigma70Definition.getSeriesAll_Unanchored(0.7));

    // Constructor
    public Parallel_v3(String referenceFile, String dir) {
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
        private List<Gene> referenceGenes;
        private NucleotideSequence nucleotides;

        // Constructor
        public GeneThreadFine(Gene gene, List<Gene> referenceGenes, NucleotideSequence nucleotides) {
            this.gene = gene;
            this.referenceGenes = referenceGenes;
            this.nucleotides = nucleotides;
        }

        // Invoke
        @Override
        public Void call() {
            for (Gene referenceGene : referenceGenes) {
                if (Homologous(gene.sequence, referenceGene.sequence)) {
                    Match prediction = PredictPromoter(GetUpstreamRegion(nucleotides, gene));
                    if (prediction != null) {
                        consensus.compute(referenceGene.name, (k,v) -> { v.addMatch(prediction); return v; });
                        consensus.compute("all", (k,v) -> {v.addMatch(prediction); return v; });
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
            // For each gene in the Ecoli file
            for (Gene gene : record.genes) {
                callableList.add(new GeneThreadFine(gene, referenceGenes, record.nucleotides));
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
