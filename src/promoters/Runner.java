package promoters;

import java.util.*;

public class Runner {
    public static void main(String[] args) throws Exception {
        //System.out.println(String.join(" ",args));
        System.out.println("PARALLEL SPEEDUP");

        // Arg parsed variables
        int maxThreads = Runtime.getRuntime().availableProcessors();
        int startThread = maxThreads;
        int iterations = 1;
        int time = 1;
        int startdelay = 0;

        // Which version do be run
        List<Parallel.Version> versions = new ArrayList<>();

        // Based off args can choose iterations to use, threads to start from
        for (int i = 0; i < args.length; i++) {

            // Add all versions
            if (args[i].equals("-all")) {
                System.out.println("-all detected");
                startThread = 0;
            }

            // Add all versions
            if (args[i].equals("-sd")) {
                startdelay = Integer.parseInt(args[++i]);
                System.out.println(String.format("-sd %d detected", startdelay));
            }

            // Add a version
            if (args[i].equals("+vall")) {
                System.out.println("+vall detected");
                versions = Arrays.asList(Parallel.Version.values());
            }
            else if (args[i].contains("+v")) {
                int version = Integer.parseInt(args[i].substring(args[i].length() - 1));
                System.out.println(String.format("+v%d detected", version));
                versions.add(Parallel.Version.values()[version - 1]);
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
            versions = Arrays.asList(Parallel.Version.values());

        // Sort
        versions = new ArrayList<>(new HashSet<>(versions));
        Collections.sort(versions);
        System.out.println("");

        // Sleep
        Thread.sleep(startdelay * 1000);

        // Setup
        String referenceFile = "../referenceGenes.list";
        String dir = "../Ecoli";
        double runtime;

        // Sequential Benchmark
        Sequential sequential = null;
        double benchmark = 0;

        // RUN IT
        System.out.println(String.format("SEQUENTIAL over %d iteration%s", iterations, (iterations > 1 ? "s" : "")));
        for (int i = 0; i < iterations; i++) {
            sequential = new Sequential(referenceFile, dir);

            // Cool the jets before spinning up again
            Thread.sleep(time * 1000);
            runtime = sequential.Run();
            benchmark +=  runtime  / iterations;
            if (iterations > 1)
                System.out.println(String.format("IT%d: %.3fs", i+1, runtime));
        }
        System.out.println(String.format("AVG: %.3fs", benchmark));
        System.out.println("");

        // Parallel Benchmark
        Parallel parallel = null;
        double runtimeAvg;

        // RUN IT
        for (Parallel.Version version : versions ) {
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
                System.out.println(String.format("AVG: %.3fs (%.3f speedup)", runtimeAvg, benchmark/runtimeAvg));
            }

            // Check that the implementation is correct
            if (!parallel.ResultsString().equals(sequential.ResultsString())) {
                sequential.DisplayResults();
                parallel.DisplayResults();
                throw new Exception("Parallelisation failed");
            }

            // Output
            //System.out.println(String.format("%s Results", version));
            //parallel.DisplayResults();
            System.out.println("");
        }
    }
}
//java -jar promoter.jar