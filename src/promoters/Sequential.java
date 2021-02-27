package promoters;

import edu.au.jacobi.pattern.Match;
import edu.au.jacobi.pattern.Series;
import jaligner.BLOSUM62;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;
import qut.*;

import java.io.*;
import java.util.*;

public class Sequential {
    private String referenceFile;
    private String dir;
    private final Matrix BLOSUM_62 = BLOSUM62.Load();
    private byte[] complement = new byte['z'];
    {
        complement['C'] = 'G'; complement['c'] = 'g'; complement['G'] = 'C'; complement['g'] = 'c';
        complement['T'] = 'A'; complement['t'] = 'a'; complement['A'] = 'T'; complement['a'] = 't';
    }

    public HashMap<String, Sigma70Consensus> consensus = new HashMap<>();
    private Series sigma70_pattern = Sigma70Definition.getSeriesAll_Unanchored(0.7);

    // Constructor
    public Sequential(String referenceFile, String dir) {
        this.referenceFile = referenceFile;
        this.dir = dir;
    }

    private List<Gene> ParseReferenceGenes(String referenceFile) throws IOException {
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

    private boolean Homologous(PeptideSequence A, PeptideSequence B) {
        return SmithWatermanGotoh.align(
                new Sequence(A.toString()),
                new Sequence(B.toString()), BLOSUM_62, 10f, 0.5f).calculateScore() >= 60;
    }

    private NucleotideSequence GetUpstreamRegion(NucleotideSequence dna, Gene gene) {
        int upStreamDistance = 250;
        if (gene.location < upStreamDistance)
            upStreamDistance = gene.location - 1;
        if (gene.strand == 1)
            return new NucleotideSequence(Arrays.copyOfRange(dna.bytes, gene.location - upStreamDistance - 1, gene.location - 1));
        else {
            byte[] result = new byte[upStreamDistance];
            int reverseStart = dna.bytes.length - gene.location + upStreamDistance;
            for (int i = 0; i < upStreamDistance; i++)
                result[i] = complement[dna.bytes[reverseStart - i]];
            return new NucleotideSequence(result);
        }
    }

    private void ProcessDir(List<String> list, File dir) {
        if (dir.exists())
            for (File file : dir.listFiles())
                if (file.isDirectory())
                    ProcessDir(list, file);
                else
                    list.add(file.getPath());
    }

    private List<String> ListGenbankFiles(String dir) {
        List<String> list = new ArrayList<String>();
        ProcessDir(list, new File(dir));
        return list;
    }

    private GenbankRecord Parse(String file) throws IOException {
        GenbankRecord record = new GenbankRecord();
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
        record.Parse(reader);
        reader.close();
        return record;
    }

    private Match PredictPromoter(NucleotideSequence upStreamRegion) {
        return BioPatterns.getBestMatch(sigma70_pattern, upStreamRegion.toString());
    }

    public void DisplayResults() {
        for (Map.Entry<String, Sigma70Consensus> entry : consensus.entrySet())
            System.out.println(entry.getKey() + " " + entry.getValue());
    }

    public String ResultsString() {
        return consensus.values().toString();
    }

    public double Run() throws IOException {
        long startTime = System.nanoTime();
        List<Gene> referenceGenes = ParseReferenceGenes(referenceFile);
        // For each Ecoli file
        for (String filename : ListGenbankFiles(dir)) {
            GenbankRecord record = Parse(filename);
            // For each gene in the reference file
            for (Gene referenceGene : referenceGenes) {
                // For each gene in the Ecoli file
                for (Gene gene : record.genes) {
                    // Nearly the entire CPU time is taken in Homologous
                    if (Homologous(gene.sequence, referenceGene.sequence)) {
                        NucleotideSequence upStreamRegion = GetUpstreamRegion(record.nucleotides, gene);
                        Match prediction = PredictPromoter(upStreamRegion);
                        if (prediction != null) {
                            // Add matches to specific and total Sigma70Consensus objects
                            consensus.get(referenceGene.name).addMatch(prediction);
                            consensus.get("all").addMatch(prediction);
                        }
                    }
                }
            }
        }
        long endTime = System.nanoTime();
        return (endTime - startTime) / 1e9;
    }

    public static void main(String[] args) throws IOException {
        // Arg parsed variables
        System.out.println("~OLD~");


        int iterations = 1;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-i")) {
                iterations = Integer.parseInt(args[++i]);
                System.out.println(String.format("-i %d detected", iterations));
            }
        }

        // Setup variables
        String referenceFile = "../referenceGenes.list";
        String dir = "../Ecoli";
        double runtime;
        Sequential sequential = null;

        System.out.print(String.format("Sequential over %d iteration%s:", iterations, (iterations > 1 ? "s" : "")));
        runtime = 0;
        for (int i = 0; i < iterations; i++) {
            sequential = new Sequential(referenceFile, dir);
            runtime += sequential.Run() / iterations;
        }
        System.out.println(String.format(" %.3fs %s", runtime, (iterations > 1 ? "avg" : "")));

        // Output
        System.out.println("Results");
        sequential.DisplayResults();
    }
}