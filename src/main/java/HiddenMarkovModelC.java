/** This class implements the topology C of Hidden Markov Model for the Ston Stemmer
 @author Matteo Lisotto (matteo.lisotto@gmail.com)
 */

import javafx.util.Pair;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Collectors;

public class HiddenMarkovModelC implements HiddenMarkovModel {
    /** Number of states */
    private int numStates;

    /** Number of prefix states */
    private int numPrefixes;

    /** Number of suffix states */
    private int numSuffixes;

    /** Length of sigma */
    private int sigmaSize;

    /** Sigma */
    private List<Character> sigma;

    /** Initial state probabilities */
    private double initialProbabilities[];

    /** Transition probabilities */
    private double transitionMatrix[][];

    /** Emission probabilities */
    private Map<Pair<Integer, Character>, Double> emissionMatrix;

    /** Smoothing variable */
    private double smooth[];

    /** Initializes an HiddenMarkovModelC.
     @param numPrefixes number of prefix states
     @param numSuffixes number of suffix states
     @param sigma the alphabet
     @param emissions emission matrix
     */
    public HiddenMarkovModelC(int numPrefixes, int numSuffixes, List<Character> sigma, Map<Pair<Integer, Character>, Double> emissions) {
        this(numPrefixes, numSuffixes, sigma);
        this.emissionMatrix = emissions;
    }

    /** Initializes an HiddenMarkovModelC.
     @param numPrefixes number of prefix states
     @param numSuffixes number of suffix states
     @param sigma the alphabet
     */
    public HiddenMarkovModelC(int numPrefixes, int numSuffixes, List<Character> sigma) {
        this.numPrefixes = numPrefixes;
        this.numSuffixes = numSuffixes;
        this.numStates = numPrefixes + numSuffixes;

        this.initialProbabilities = new double[numStates];
        this.transitionMatrix = new double[numStates][numStates];
        this.emissionMatrix = new LinkedHashMap<>();
        this.sigma = new ArrayList<>(sigma);

        /* Initialize the parameter initialProbability */
        double initialProbability = 1.0 / (numPrefixes - 2);
        Arrays.fill(this.initialProbabilities, 0, numPrefixes - 2, initialProbability);

        for (int i = 0; i < numPrefixes - 2; i++) {
            this.transitionMatrix[i][i] = 0.5;
            this.transitionMatrix[i][i+1] = 0.5;
        }

        for (int i = numPrefixes - 2; i < numStates - 1; i++) {
            this.transitionMatrix[i][i + 1] = 1;
        }

        /* Initializes the emission probabilities */
        this.sigmaSize = sigma.size();
        double emissionProbability = 1.0 / sigmaSize;


        for (int i = 0; i < numStates; i++) {
            for(Character c: sigma) {
                emissionMatrix.put(new Pair<>(i, c), emissionProbability);
            }
        }

    }

    /** Trains the Hidden Markov Model using the Baum-Welch Algorithms.
     @param observations the training set
     */
    public <T> void train (T observations) {
        train(observations, 1);
    }
    /** Trains the Hidden Markov Model using the Baum-Welch Algorithms.
     @param observations the training set
     @param steps the number of steps
     */
    public <T> void train(T observations, int steps) {
        double forward [][];
        double backward [][];

        double initialP [] = new double[numStates];
        double transitionM [][] = new double[numStates][numStates];
        Map<Pair<Integer, Character>, Double> emissionM = new LinkedHashMap<>();

        // WIP
        String observation = ArrayList.class.cast(observations).stream().map(Object::toString)
                .collect(Collectors.joining("")).toString();
                //String.join("", ArrayList.class.cast(observations));


        for (int s = 0; s < steps; s++) {
        //for (Object obs: ArrayList.class.cast(observations)) {
                //String observation = obs.toString();

            /* Calculation of Forward and Backward variables */
                forward = forwardProc(observation);
                backward = backwardProc(observation);

            /* Re-estimation of initial state probabilities */
                for (int i = 0; i < numStates; i++)
                    initialP[i] = gamma(i, 0, forward, backward);

            /* Re-estimation of transition probabilities */
                for (int i = 0; i < numStates; i++) {
                    for (int j = 0; j < numStates; j++) {
                        transitionM[i][j] = epsilon(i, j, observation, forward, backward);
                    }
                }
      
            /* Re-estimation of emission probabilities */
                for (int i = 0; i < numStates; i++) {
                    for (int k = 0; k < sigmaSize; k++) {
                        double num = 0;
                        double denom = 0;
                        for (int t = 0; t <= observation.length() - 1; t++) {
                            double g = gamma(i, t, forward, backward);
                            num += g * (sigma.get(k).equals(observation.charAt(t)) ? 1 : 0);
                            denom += g;
                        }
                        emissionM.put(new Pair<>(i, sigma.get(k)), divide(num, denom));
                    }
                }
                initialProbabilities = initialP;
                transitionMatrix = transitionM;
                emissionMatrix.putAll(emissionM);
            }
    }


    /** Calculation of Forward-Variables f(i,t) for state i at time
     t for output sequence O with the current HiddenMarkovModelC parameters
     @param observation the output sequence O
     @return an array f(i,t) over states and times, containing
     the Forward-variables.
     */
    private double[][] forwardProc(String observation) {
        int lenObs = observation.length();
        double[][] forward = new double[numStates][lenObs];
        smooth = new double[lenObs];

    /* Initialization (time 0) */
        for (int i = 0; i < numStates; i++) {
            forward[i][0] = initialProbabilities[i] * emissionMatrix.get(new Pair<>(i, observation.charAt(0)));
            smooth[0] += forward[i][0];
        }

        for (int i = 0; i < numStates; i++)
            forward[i][0] = divide(forward[i][0], smooth[0]);

    /* Induction */
        for (int t = 0; t <= lenObs - 2; t++) {
            for (int j = 0; j < numStates; j++) {
                forward[j][t+1] = 0;
                for (int i = 0; i < numStates; i++)
                    forward[j][t+1] += (forward[i][t] * transitionMatrix[i][j]);
                forward[j][t+1] *= emissionMatrix.get(new Pair<>(j, observation.charAt(t + 1)));
                smooth[t+1] += forward[j][t+1];
            }

            for(int i = 0; i < numStates; i++)
                forward[i][t+1] = divide(forward[i][t+1], smooth[t+1]);
        }

        return forward;
    }

    /** Calculation of  Backward-Variables b(i,t) for state i at time
     t for output sequence O with the current HiddenMarkovModelC parameters
     @param observable the output sequence O
     @return an array b(i,t) over states and times, containing
     the Backward-Variables.
     */
    private double[][] backwardProc(String observable) {
        int lenObs = observable.length();
        double backward[][] = new double[numStates][lenObs];
        
    /* Initialization (time 0) */
        for (int i = 0; i < numStates; i++)
            backward[i][lenObs - 1] = 1;

    /* Induction */
        for (int t = lenObs - 2; t >= 0; t--) {
            for (int i = 0; i < numStates; i++) {
                backward[i][t] = 0;
                for (int j = 0; j < numStates; j++)
                    backward[i][t] += (backward[j][t+1] * transitionMatrix[i][j] * emissionMatrix.get(new Pair<>(j, observable.charAt(t + 1))));
                backward[i][t] = divide(backward[i][t], smooth[t+1]);
            }
        }

        return backward;
    }

    /** Calculation of probability P(X_t = s_i, X_t+1 = s_j | O, m).
     @param i the number of state s_i
     @param j the number of state s_j
     @param observation an output sequence o
     @param forward the Forward-Variables
     @param backward the Backward-Variables
     @return epsilon of the state s_i
     */
    private double epsilon(int i, int j, String observation, double[][] forward, double[][] backward) {
        double num = 0;
        double denom = 0;
        int lenObs = observation.length();

        for(int t = 0; t < lenObs - 1; t++) {
            num += forward[i][t] * transitionMatrix[i][j] * emissionMatrix.get(new Pair<>(j, observation.charAt(t + 1))) * backward[j][t+1];
        }

        for(int t = 0; t < lenObs - 1; t++) {
            for(int k = 0; k < numStates; k++) {
                denom += forward[i][t] * transitionMatrix[i][k] * emissionMatrix.get(new Pair<>(k, observation.charAt(t + 1))) * backward[k][t+1];
            }
        }
        return divide(num, denom);
    }

    /** Computes gamma(i, t)
     @param i the number of state s_i
     @param t time t
     @param forward the Forward-Variable
     @param backward the Backward-Variable
     @return gamma of the state s_i
     */
    private double gamma(int i, int t, double[][] forward, double[][] backward) {
        double num = forward[i][t] * backward[i][t];
        double denom = 0;

        for (int j = 0; j < numStates; j++)
            denom += forward[j][t] * backward[j][t];

        return divide(num, denom);
    }

    /** Prints all the parameters of an HiddenMarkovModelC */
    public void print() {
        DecimalFormat fmt = new DecimalFormat();
        fmt.setMinimumFractionDigits(20);
        fmt.setMaximumFractionDigits(20);

        System.out.println("Initial Probabilities");
        for (int i = 0; i < numStates; i++)
            System.out.println("State(" + i + ") = " + fmt.format(initialProbabilities[i]));
        System.out.println();

        System.out.println("Transition Matrix");
        for (int i = 0; i < numStates; i++) {
            for (int j = 0; j < numStates; j++)
                System.out.print("State(" + i + "," + j + ") = " +
                        fmt.format(transitionMatrix[i][j]) + "  ");
            System.out.println();
        }

        System.out.println();
        System.out.println("Emission Matrix");
        emissionMatrix.forEach((k,v) -> System.out.println("Key: " + k + " Value: " + v));
    }

    /** divides two doubles. 0 / 0 = 0! */
    private double divide(double n, double d) {
        if (n == 0)
            return 0;
        else
            return n / d;
    }
}
