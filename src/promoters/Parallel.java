package promoters;

import edu.au.jacobi.pattern.Match;
import jaligner.BLOSUM62;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.matrix.Matrix;
import qut.*;

import java.io.*;
import java.util.*;

public abstract  class Parallel {
    protected String referenceFile;
    protected String dir;
    protected final Matrix BLOSUM_62 = BLOSUM62.Load();
    protected byte[] complement = new byte['z'];
    {
        complement['C'] = 'G'; complement['c'] = 'g'; complement['G'] = 'C'; complement['g'] = 'c';
        complement['T'] = 'A'; complement['t'] = 'a'; complement['A'] = 'T'; complement['a'] = 't';
    }

    // Constructor
    public Parallel(String referenceFile, String dir) {
        this.referenceFile = referenceFile;
        this.dir = dir;
    }

    protected void ProcessDir(List<String> list, File dir) {
        if (dir.exists())
            for (File file : dir.listFiles())
                if (file.isDirectory())
                    ProcessDir(list, file);
                else
                    list.add(file.getPath());
    }

    protected boolean Homologous(PeptideSequence A, PeptideSequence B) {
        return SmithWatermanGotoh.align(new Sequence(A.toString()), new Sequence(B.toString()), BLOSUM_62, 10f, 0.5f).calculateScore() >= 60;
    }

    protected NucleotideSequence GetUpstreamRegion(NucleotideSequence dna, Gene gene) {
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

    protected List<String> ListGenbankFiles(String dir) {
        List<String> list = new ArrayList<String>();
        ProcessDir(list, new File(dir));
        return list;
    }

    protected GenbankRecord Parse(String file) throws IOException {
        GenbankRecord record = new GenbankRecord();
        BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
        record.Parse(reader);
        reader.close();
        return record;
    }


    // ABSTRACT METHODS
    protected abstract List<Gene> ParseReferenceGenes(String referenceFile) throws IOException;

    public abstract void DisplayResults();

    protected abstract Match PredictPromoter(NucleotideSequence upStreamRegion);

    protected abstract void Run(Integer threads) throws IOException;

    public abstract String ResultsString();


    // MAIN USAGE

    public double Time(Integer threads) throws IOException {
        long startTime = System.nanoTime();

        // Abstract method
        Run(threads);

        long endTime = System.nanoTime();
        return (endTime - startTime) / 1e9;
    }

    public enum Version { VERSION1, VERSION2, VERSION3 }

    public static void main(String[] args) throws IOException, InterruptedException {
        System.out.println("PARALLEL COMPARISON");

        // Arg parsed variables
        int maxThreads = Runtime.getRuntime().availableProcessors();
        int startThread = maxThreads;
        int iterations = 1;
        int time = 1;

        // Which version do be run
        List<Version> versions = new ArrayList<>();

        // Based off args can choose iterations to use, threads to start from
        for (int i = 0; i < args.length; i++) {
            // Add all versions
            if (args[i].equals("-all")) {
                System.out.println("-all detected");
                startThread = 0;
            }

            // Add a version
            if (args[i].equals("+vall")) {
                System.out.println("+vall detected");
                versions = Arrays.asList(Version.values());
            }
            else if (args[i].contains("+v")) {
                int version = Integer.parseInt(args[i].substring(args[i].length() - 1));
                System.out.println(String.format("+v%d detected", version));
                versions.add(Version.values()[version - 1]);
            }

            // How many iteratiosn
            if (args[i].equals("-i")) {
                iterations = Integer.parseInt(args[++i]);
                System.out.println(String.format("-i %d detected", iterations));
            }

            // How long to wait between runs
            if (args[i].equals("-t")) {
                time = Integer.parseInt(args[++i]);
                System.out.println(String.format("-t %d detected", time));
            }
        }

        // Make sure args had something
        if (versions.isEmpty())
            versions = Arrays.asList(Version.values());

        // Sort
        versions = new ArrayList<>(new HashSet<>(versions));
        Collections.sort(versions);
        System.out.println();

        // Setup variables
        String referenceFile = "../referenceGenes.list";
        String dir = "../Ecoli";
        double runtime;
        double runtimeAvg;
        Parallel parallel = null;


        // RUN IT
        for (Version version :versions ) {
            for (int threads = startThread; threads <= maxThreads; threads += 2) {

                // For a verion and a thread, average runtime over iterations
                System.out.println(String.format("%s using %d threads over %d iteration%s",
                        version, (threads < 1 ? 1 : threads), iterations, (iterations > 1 ? "s" : "")));

                runtimeAvg = 0;
                for (int i = 0; i < iterations; i++) {
                    switch (version) {
                        case VERSION1:
                            parallel = new Parallel_v1(referenceFile, dir);
                            break;
                        case VERSION2:
                            parallel = new Parallel_v2(referenceFile, dir);
                            break;
                        case VERSION3:
                            parallel = new Parallel_v3(referenceFile, dir);
                            break;

                        default:
                            throw new IllegalStateException("Unexpected value: " + version);
                    }

                    // Cool the jets before spinning up again
                    Thread.sleep(time * 1000);
                    runtime = parallel.Time(threads < 1 ? 1 : threads);
                    runtimeAvg += runtime / iterations;
                    if (iterations > 1)
                        System.out.println(String.format("IT%d: %.3fs", i+1, runtime));
                }

                System.out.println(String.format("AVG: %.3fs", runtimeAvg));
            }

            // Output
            System.out.println(String.format("%s Results", version));
            parallel.DisplayResults();
            System.out.println("");
        }
    }
}
