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
import java.util.concurrent.*;

public class Parallel_v1 extends Parallel {

    // Threading
    protected HashMap<String, Sigma70Consensus> consensus = new HashMap<>();
    public Series sigma70_pattern = Sigma70Definition.getSeriesAll_Unanchored(0.7);

    // Constructor
    public Parallel_v1(String referenceFile, String dir) {
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
        return BioPatterns.getBestMatch(sigma70_pattern, upStreamRegion.toString());
    }

    private class GenesHomologous {
        public Gene gene;
        public Gene referenceGene;
        public boolean homologous;
        private NucleotideSequence nucleotides;

        //Constructor
        public GenesHomologous(Gene gene, Gene referenceGene, boolean homologous, NucleotideSequence nucleotideSequence) {
            this.gene = gene;
            this.referenceGene = referenceGene;
            this.homologous = homologous;
            this.nucleotides = nucleotideSequence;
        }
    }

    private class GeneThreadFine implements Callable<GenesHomologous> {
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
        public GenesHomologous call() {
            return new GenesHomologous(gene, referenceGene,
                    Homologous(gene.sequence, referenceGene.sequence), nucleotides);
        }
    }

    @Override
    public void Run(Integer threads) throws IOException {
        // Scheduling setup
        ExecutorService executorService = Executors.newFixedThreadPool(threads);
        List<Callable<GenesHomologous>> callableList = new ArrayList<>();

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
            List<Future<GenesHomologous>> resultList = executorService.invokeAll(callableList);
            for (int i = 0; i < resultList.size(); i++) {
                Future<GenesHomologous> future = resultList.get(i);
                GenesHomologous result = future.get();
                if (result.homologous) {
                    NucleotideSequence upStreamRegion = GetUpstreamRegion(result.nucleotides, result.gene);
                    Match prediction = PredictPromoter(upStreamRegion);
                    if (prediction != null) {
                        // Add matches to specific and total Sigma70Consensus objects
                        consensus.get(result.referenceGene.name).addMatch(prediction);
                        consensus.get("all").addMatch(prediction);
                    }
                }
            }
            executorService.shutdown();
        } catch (InterruptedException| ExecutionException ex) {
            ex.printStackTrace();
        }
    }


}
