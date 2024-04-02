package com.company;

import com.mathworks.engine.*;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutionException;

import static com.company.SequenceClass.maxsearch;

class SequenceClass implements Comparable<SequenceClass> {
    // Holds sequence information in StringBuilder form
    public StringBuilder sequenceBuilder;
    // Holds sequence information in String form
    public String sequenceString;
    // Holds sequence "magnitude"
    public double magnitude;
    // Searching for maximum (change to -1 if not)
    public static int maxsearch = 1;

    /**
     * Constructor for sequence represented by StringBuilder
     *
     * @param sequence Sequence in StringBuilder format
     * @param desiredIndex Desired index for processing sequence
     * @param processingMethod Method to process sequence
     * @param eng MatlabEngine to execute method
     */
    public SequenceClass(StringBuilder sequence, int desiredIndex, String processingMethod, MatlabEngine eng) throws ExecutionException, InterruptedException {
        sequenceBuilder = new StringBuilder(sequence);
        magnitude = processSequence(sequenceBuilder.toString(), desiredIndex, processingMethod, eng);
    }

    /**
     * Constructor for sequence represented by String
     *
     * @param sequence Sequence in String format
     * @param desiredIndex Desired index for processing sequence
     * @param processingMethod Method to process sequence
     * @param eng MatlabEngine to execute method
     */
    public SequenceClass(String sequence, int desiredIndex, String processingMethod, MatlabEngine eng) throws ExecutionException, InterruptedException {
        sequenceString = sequence;
        magnitude = processSequence(sequenceString, desiredIndex, processingMethod, eng);
    }

    /**
     * Process the sequence string using MatlabEngine
     *
     * @param sequence String representation of the sequence
     * @param desiredIndex Desired index for processing sequence
     * @param processingMethod Method to process sequence
     * @param eng MatlabEngine to execute method
     * @return Sequence magnitude
     */
    private double processSequence(String sequence, int desiredIndex, String processingMethod, MatlabEngine eng) throws ExecutionException, InterruptedException {
        if (sequence.length() > 1) {
            return eng.feval(processingMethod, sequence.toCharArray(), desiredIndex);
        } else {
            return 0.0;
        }
    }

    /**
     * Compare method for Sequence class
     * @param otherSequence Sequence to compare with
     * @return Comparison result
     */
    public int compareTo(SequenceClass otherSequence) {
        if (maxsearch==1)
        {
            return Double.compare(magnitude, otherSequence.magnitude);
        }
        else
        {
            return Double.compare(otherSequence.magnitude, magnitude);
        }
    }
}

public class MainSearch {
    // Global counter for sequence rejections
    public static int rejectionCounter = 0;
    // Instance of MatlabEngine
    public static MatlabEngine eng;
    public static String modeltype="runningforjava"; // Change to runningforjavacgDNA for the cgDNA parametrization (default is crystal
    
    public static void main(String[] args) throws InterruptedException, ExecutionException, IOException {
        eng = MatlabEngine.startMatlab(); // Start Matlab engine
        eng.eval("addpath('" + System.getProperty("user.dir") + "')"); // Add working directory to Matlab engine's path
        String[] metrics= {
            "angoffset",
                    "persistence_length",
                    "angdev",
                    "bend",
                    "extension",
                    "contour_length",
                    "beam_projection",
                    "curvsum",
                    "overall_compliance",
                    "bending_compliance"
        };
        // Loop for each desired index from 1 to 10, each indexing possible metric
        for (int desiredIndex = 1; desiredIndex <= 9; desiredIndex++) {
            ArrayList<SequenceClass> DNASequences = new ArrayList<>();
            ArrayList<String> bases = new ArrayList<>();
            double baseRate;
            LinkedList<Double> enterOrder = new LinkedList<>();
            ArrayList<Double> rankedOrder = new ArrayList<>();

            // Initialize rankedOrder and enterOrder with zeros
            for (int k = 0; k < 100; k++) {
                rankedOrder.add(0.0);
                enterOrder.add(0.0);
            }

            // Add bases
            bases.add("A");
            bases.add("C");
            bases.add("G");
            bases.add("T");

            ArrayList<SequenceClass> builder = new ArrayList<>();
            builder.add(new SequenceClass("", desiredIndex, modeltype, eng));

            // Start building sequences
            for (int d = 0; d < 100; d++) {
                int sup = builder.size();
                System.out.println(d + "/" + 100);
                for (int c = 0; c < sup; c++) {
                    extendFirstSequence(builder, desiredIndex, modeltype);
                }
                Collections.sort(builder);
                if (d > 0) {
                    removeWeakestSequences(builder);
                }
            }

            StringBuilder startingSequence = new StringBuilder(builder.get(builder.size() - 1).sequenceString);
            DNASequences.add(new SequenceClass(startingSequence, desiredIndex, modeltype, eng));

            // Setup variables for acceptance-rejection test
            int runs = 10;
            double changedBasesGrowth;
            double baseStep;
            double tempRate = 3.49651;
            double tempStep = 1.10866 / (runs * runs);
            double maxChangedBases;
            changedBasesGrowth = 0.5;
            maxChangedBases = 33;
            baseRate = Math.pow(maxChangedBases - 1, (1 / changedBasesGrowth)) / runs;
            baseStep = baseRate;
            double acceptanceParam;
            double rejectionNumber = 0;
            double acceptanceNumber = 0;
            int counter = 0;

            // Main loop for changing bases
            for (double currentBaseChangeAmount = maxChangedBases; currentBaseChangeAmount > 1; currentBaseChangeAmount = maxChangedBases - Math.pow(baseRate, changedBasesGrowth)) {
                baseRate = baseRate + baseStep;
                for (int innerBaseIteration = 0; innerBaseIteration < runs; innerBaseIteration++) {
                    acceptanceParam = changeBase(DNASequences, bases, rankedOrder, enterOrder, (int) Math.floor(Math.exp(tempRate)), currentBaseChangeAmount, desiredIndex);
                    if (acceptanceParam > 0)
                        acceptanceNumber++;
                    tempRate = tempRate + tempStep;

                }
                System.out.println(counter);
                counter = counter + 1;
            }

            Collections.sort(DNASequences);
            // Writing the top sequence to a file
            for (int b = 1; b < 2; b++) {
                FileWriter myWriter = new FileWriter("SampleMax.txt", true);
                myWriter.write(metrics[desiredIndex-1]+'\n');
                myWriter.write(DNASequences.get(DNASequences.size() - b).magnitude + " " + DNASequences.get(DNASequences.size() - b).sequenceBuilder + '\n');
                myWriter.close();
            }
        }
        eng.close(); // Close Matlab engine
    }
    /**
     * Mutates the last Sequence in the DNASequences list by changing the bases at random positions.
     * The number of bases changed is determined by baseChangeNumber. The method returns the number
     * of rejections encountered before the change was accepted, and it updates the ranked and keytime lists.
     *
     * @param dnaSequences List of Sequence instances
     * @param baseOptions List of possible bases: 'A', 'C', 'G', 'T'
     * @param ranked List of magnitudes of previously accepted changes, in descending order
     * @param keytime List of magnitudes of the most recent changes, in the order they were accepted
     * @param percentile Index to decide the acceptance of the change
     * @param baseChangeNumber Average number of bases to change in the mutation
     * @param desiredIndex Index used by the Sequence constructor
     * @return Number of rejections before the change was accepted
     * @throws ExecutionException, InterruptedException
     */
    public static double changeBase (ArrayList<SequenceClass> dnaSequences, ArrayList<String> baseOptions, ArrayList<Double> ranked, LinkedList<Double> keytime, int percentile, double baseChangeNumber, int desiredIndex) throws ExecutionException, InterruptedException {

        int rejectionCountBeforeAcceptance;
        // Round baseChangeNumber to the nearest integer, with probabilistic rounding
        if (Math.random() > (baseChangeNumber - Math.floor(baseChangeNumber)))
            baseChangeNumber = Math.floor(baseChangeNumber);
        else
            baseChangeNumber = Math.ceil(baseChangeNumber);

        StringBuilder currentSequence = new StringBuilder(dnaSequences.get(dnaSequences.size() - 1).sequenceBuilder);

        // Change bases at random positions in the Sequence
        for (int basesToChangeCount = 0; basesToChangeCount < baseChangeNumber; basesToChangeCount++) {
            int positionToChange = (int) (Math.random() * dnaSequences.get(0).sequenceBuilder.length());
            String baseAtPosition = Character.toString(dnaSequences.get(dnaSequences.size() - 1).sequenceBuilder.charAt(positionToChange));
            int baseIndex = (int) (Math.random() * 3);

            baseOptions.remove(baseAtPosition);
            currentSequence.replace(positionToChange, positionToChange + 1, baseOptions.get(baseIndex));
            baseOptions.add(baseAtPosition);
        }

        SequenceClass mutatedSequence = new SequenceClass(currentSequence, desiredIndex, modeltype, eng);

        // Accept the change if it improves the Sequence or is better than the percentile in ranked
        if (mutatedSequence.magnitude*maxsearch > dnaSequences.get(dnaSequences.size() - 1).magnitude*maxsearch || mutatedSequence.magnitude*maxsearch > ranked.get(percentile)*maxsearch) {
            dnaSequences.add(mutatedSequence);
            rejectionCountBeforeAcceptance = rejectionCounter;
            rejectionCounter = 0;
            updateDistribution(ranked,keytime,mutatedSequence.magnitude);
            return rejectionCountBeforeAcceptance;
        } else {
            rejectionCounter++;
            updateDistribution(ranked,keytime,mutatedSequence.magnitude);
            return 0;
        }
    }


    /**
     * Update ranked and keytime lists with the new value.
     * The ranked list is sorted in ascending order, with the highest value at the beginning of the list.
     * The keytime list acts as a queue of the most recent changes.
     *
     * @param ranked List of magnitudes of previously accepted changes, in descending order
     * @param keytime List of magnitudes of the most recent changes, in the order they were accepted
     * @param newValue The magnitude of the new change that needs to be added to both lists
     */
    public static void updateDistribution(ArrayList<Double> ranked, LinkedList<Double> keytime, double newValue) {
        // Remove the oldest change from keytime and ranked
        double oldestChange = keytime.removeLast();
        keytime.addFirst(newValue);
        ranked.remove(oldestChange);

        // Add the new value to the correct position in ranked, keeping the list sorted
        for (int index = 0; index < ranked.size(); index++) {
            if (ranked.get(index)*maxsearch > newValue*maxsearch) {
                ranked.add(index, newValue);
                return;
            }
        }

        // If the new value is the smallest, add it at the end
        ranked.add(newValue);
    }
    /**
     * Removes the weakest 75% of the sequences from the given list.
     * The input list should be sorted in ascending order of 'strength' or 'fitness'.
     *
     * @param sequences A list of Sequence objects sorted in ascending order of their fitness
     */
    public static void removeWeakestSequences(ArrayList<SequenceClass> sequences) {
        int sequencesToRemove = (int)(sequences.size() * 0.75);

        // Remove the weakest 75% sequences from the list.
        for (int index = 0; index < sequencesToRemove; index++) {
            sequences.remove(0);
        }
    }

    /**
     * Extend the first Sequence in the list with each of the four bases (A, C, G, T),
     * then add each extended Sequence to the list.
     * Finally, remove the original (non-extended) Sequence from the list.
     *
     * @param sequenceList A list of Sequence objects
     * @param chosenIndex The index to be set for the new Sequences
     * @param directive The directive to be set for the new Sequences
     */
    public static void extendFirstSequence(ArrayList<SequenceClass> sequenceList, int chosenIndex, String directive) throws ExecutionException, InterruptedException {
        // Get the code of the first Sequence in the list
        String originalCode = sequenceList.get(0).sequenceString;

        // Extend the original code with each base and add each extended Sequence to the list
        sequenceList.add(new SequenceClass(originalCode + "A", chosenIndex, directive, eng));
        sequenceList.add(new SequenceClass(originalCode + "C", chosenIndex, directive, eng));
        sequenceList.add(new SequenceClass(originalCode + "G", chosenIndex, directive, eng));
        sequenceList.add(new SequenceClass(originalCode + "T", chosenIndex, directive, eng));

        // Remove the original Sequence from the list
        sequenceList.remove(0);
    }

}
