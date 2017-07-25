/**
 * This class implements the topology C of Hidden Markov Model for the Ston Stemmer
 *
 * @author Matteo Lisotto (matteo.lisotto@gmail.com)
 */
import java.util.Arrays;
import java.util.Map;
import java.util.LinkedHashMap;
import javafx.util.Pair;
import java.text.DecimalFormat;
import java.util.stream.Collectors;
import java.util.List;
import java.util.ArrayList;

public class HiddenMarkovModelCSmoothed implements HiddenMarkovModel {
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
  private List<Double[]> smooth;

  /** Emissions sums */
  private double emissionsSums[];

  /**
   * Initializes an HiddenMarkovModelCSmoothed.
   *
   * @param numPrefixes number of prefix states
   * @param numSuffixes number of suffix states
   * @param sigma the alphabet
   */
  public HiddenMarkovModelCSmoothed(int numPrefixes, int numSuffixes, List<Character> sigma) {
    this.numPrefixes = numPrefixes;
    this.numSuffixes = numSuffixes;
    this.numStates = numPrefixes + numSuffixes;

    this.initialProbabilities = new double[numStates];
    this.transitionMatrix = new double[numStates][numStates];
    this.emissionMatrix = new LinkedHashMap<>();
    this.sigma = new ArrayList<>(sigma);
    this.smooth = new ArrayList<Double[]>();

    /* Initialize the parameter initialProbability */
    double initialProbability = 1.0 / (numPrefixes - 2);
    Arrays.fill(this.initialProbabilities, 0, numPrefixes - 2, initialProbability);

    for (int i = 0; i < numPrefixes - 2; i++) {
      this.transitionMatrix[i][i] = 0.5;
      this.transitionMatrix[i][i + 1] = 0.5;
    }

    for (int i = numPrefixes - 2; i < numStates - 1; i++) {
      this.transitionMatrix[i][i + 1] = 1;
    }

    /* Initializes the emission probabilities */
    this.sigmaSize = sigma.size();
    double emissionProbability = 1.0 / sigmaSize;

    for (int i = 0; i < numStates; i++) {
      for (Character c : sigma)
        emissionMatrix.put(new Pair<Integer, Character>(i, c), emissionProbability);
    }
  }

  /**
   * Trains the Hidden Markov Model using the Baum-Welch Algorithms.
   *
   * @param observations the training set
   */
  public void train(List<String> observations) {
    train(observations, 1);
  }
  /**
   * Trains the Hidden Markov Model using the Baum-Welch Algorithms.
   *
   * @param observations the training set
   * @param steps the number of steps
   */
  public void train(List<String> observations, int steps) {
    List<Double[][]> forward;
    List<Double[][]> backward;
    int numObservations = observations.size();

    double initialP[] = new double[numStates];
    double transitionM[][] = new double[numStates][numStates];
    Map<Pair<Integer, Character>, Double> emissionM = new LinkedHashMap<>();

    for (int s = 0; s < steps; s++) {
      /* Calculation of Forward and Backward variables */
      forward = forwardProc(observations);
      backward = backwardProc(observations);

      /* Re-estimation of initial state probabilities */
      for (int i = 0; i < numStates; i++) {
        for (int k = 0; k < numObservations; k++)
          initialP[i] += gamma(i, 0, forward.get(k), backward.get(k));
        initialP[i] = divide(initialP[i], numObservations);
      }

      /* Re-estimation of transition probabilities */
      for (int i = 0; i < numStates; i++) {
        for (int j = 0; j < numStates; j++) {
          double num = 0;
          double den = 0;
          for (int k = 0; k < numObservations; k++) {
            String observation = observations.get(k);
            int lenObs = observation.length();
            for (int t = 0; t < lenObs - 1; t++) {
              num +=
                  forward.get(k)[i][t]
                      * transitionMatrix[i][j]
                      * emissionMatrix.get(
                          new Pair<Integer, Character>(j, observation.charAt(t + 1)))
                      * backward.get(k)[j][t + 1];
              for (int z = 0; z < numStates; z++)
                den +=
                    forward.get(k)[i][t]
                        * transitionMatrix[i][z]
                        * emissionMatrix.get(
                            new Pair<Integer, Character>(z, observation.charAt(t + 1)))
                        * backward.get(k)[z][t + 1];
            }
          }
          transitionM[i][j] = divide(num, den);
        }
      }

      /* Re-estimation of emission probabilities */
      for (int i = 0; i < numStates; i++) {
        for (int k = 0; k < sigmaSize; k++) {
          double num = 0;
          double denom = 0;
          for (int j = 0; j < numObservations; j++) {
            String observation = observations.get(j);
            for (int t = 0; t <= observation.length() - 1; t++) {
              double g = gamma(i, t, forward.get(j), backward.get(j));
              num += g * (sigma.get(k).equals(observation.charAt(t)) ? 1 : 0);
              denom += g;
            }
          }
          emissionM.put(new Pair<Integer, Character>(i, sigma.get(k)), divide(num, denom));
        }
      }
      initialProbabilities = initialP;
      transitionMatrix = transitionM;
      emissionMatrix.putAll(emissionM);
    }
  }

  /**
   * Decodes the Trains the Hidden Markov Model using the Baum-Welch Algorithms.
   *
   * @param observation the string to stem
   * @return string the stemmed string
   */
  public String decode(String observation) {
    int lenObs = observation.length();
    double delta[][] = new double[numStates][lenObs];
    int path[] = new int[lenObs];
    double maxProb = -100000;
    int stemPosition;

    Arrays.fill(path, -1);

    /* Calculate delta at time 0. At first finds the maxProb */
    delta[0][0] =
        Math.log(initialProbabilities[0])
            + Math.log(emissionMatrix.get(new Pair<Integer, Character>(0, observation.charAt(0))));
    maxProb = delta[0][0];
    path[0] = 0;
    for (int i = 1; i < numStates; i++) {
      delta[i][0] =
          Math.log(initialProbabilities[i])
              + Math.log(
                  emissionMatrix.get(new Pair<Integer, Character>(i, observation.charAt(0))));
      if (delta[i][0] > maxProb) {
        path[0] = i;
        maxProb = delta[i][0];
      }
    }

    /* Induction */
    for (int t = 1; t < lenObs; t++) {
      delta[0][t] =
          delta[0][t - 1]
              + Math.log(transitionMatrix[0][0])
              + Math.log(
                  emissionMatrix.get(new Pair<Integer, Character>(0, observation.charAt(t))));
      path[t] = 0;
      for (int i = 1; i < numStates; i++) {
        double prob =
            delta[i][t - 1]
                + Math.log(transitionMatrix[i][0])
                + Math.log(
                    emissionMatrix.get(new Pair<Integer, Character>(0, observation.charAt(t))));
        if (prob > delta[0][t]) {
          delta[0][t] = prob;
        }
      }
      maxProb = delta[0][t];
      path[t] = 0;

      for (int j = 1; j < numStates; j++) {
        delta[j][t] =
            delta[0][t - 1]
                + Math.log(transitionMatrix[0][j])
                + Math.log(
                    emissionMatrix.get(new Pair<Integer, Character>(j, observation.charAt(t))));
        for (int i = 1; i < numStates; i++) {
          double prob =
              delta[i][t - 1]
                  + Math.log(transitionMatrix[i][j])
                  + Math.log(
                      emissionMatrix.get(new Pair<Integer, Character>(j, observation.charAt(t))));
          if (prob > delta[j][t]) {
            delta[j][t] = prob;
          }
        }
        if (delta[j][t] > maxProb) path[t] = j;
      }
    }

    stemPosition = Arrays.asList(path).indexOf(numPrefixes);
    return stemPosition != -1 ? observation.substring(0, numPrefixes) : observation;
  }

  /**
   * Calculation of Forward-Variables f(i,t) for state i at time t for output sequence O with the
   * current HiddenMarkovModelCSmoothed parameters
   *
   * @param observation the output sequence O
   * @return an array f(i,t) over states and times, containing the Forward-variables.
   */
  private List<Double[][]> forwardProc(List<String> observations) {
    List<Double[][]> forward = new ArrayList<Double[][]>();

    for (String observation : observations) {
      int lenObs = observation.length();
      Double forwardI[][] = new Double[numStates][lenObs];
      Double smoothI[] = new Double[lenObs];

      /* Initialization (time 0) */
      for (int i = 0; i < numStates; i++) {
        forwardI[i][0] =
            initialProbabilities[i]
                * emissionMatrix.get(new Pair<Integer, Character>(i, observation.charAt(0)));
        smoothI[0] += forwardI[i][0];
      }
      for (int i = 0; i < numStates; i++) forwardI[i][0] = divide(forwardI[i][0], smoothI[0]);

      /* Induction */
      for (int t = 0; t <= lenObs - 2; t++) {
        for (int j = 0; j < numStates; j++) {
          forwardI[j][t + 1] = 0.0;
          for (int i = 0; i < numStates; i++)
            forwardI[j][t + 1] += (forwardI[i][t] * transitionMatrix[i][j]);
          forwardI[j][t + 1] *=
              emissionMatrix.get(new Pair<Integer, Character>(j, observation.charAt(t + 1)));
          smoothI[t + 1] += forwardI[j][t + 1];
        }
        for (int i = 0; i < numStates; i++)
          forwardI[i][t + 1] = divide(forwardI[i][t + 1], smoothI[t + 1]);
      }
      forward.add(forwardI);
      smooth.add(smoothI);
    }
    return forward;
  }

  /**
   * Calculation of Backward-Variables b(i,t) for state i at time t for output sequence O with the
   * current HiddenMarkovModelCSmoothed parameters
   *
   * @param observation the output sequence O
   * @return an array b(i,t) over states and times, containing the Backward-Variables.
   */
  private List<Double[][]> backwardProc(List<String> observations) {
    List<Double[][]> backward = new ArrayList<Double[][]>();

    for (int k = 0; k < observations.size(); k++) {
      String observation = observations.get(k);
      int lenObs = observation.length();
      Double backwardI[][] = new Double[numStates][lenObs];

      /* Initialization (time 0) */
      for (int i = 0; i < numStates; i++) backwardI[i][lenObs - 1] = 1.0;

      /* Induction */
      for (int t = lenObs - 2; t >= 0; t--) {
        for (int i = 0; i < numStates; i++) {
          backwardI[i][t] = 0.0;
          for (int j = 0; j < numStates; j++)
            backwardI[i][t] +=
                (backwardI[j][t + 1]
                    * transitionMatrix[i][j]
                    * emissionMatrix.get(
                        new Pair<Integer, Character>(j, observation.charAt(t + 1))));
          backwardI[i][t] = divide(backwardI[i][t], smooth.get(k)[t + 1]);
        }
      }
      backward.add(backwardI);
    }

    return backward;
  }

  /**
   * Computes gamma(i, t)
   *
   * @param i the number of state s_i
   * @param t time t
   * @param forward the Forward-Variable
   * @param backward the Backward-Variable
   * @return gamma of the state s_i
   */
  private double gamma(int i, int t, Double[][] forward, Double[][] backward) {
    double num = forward[i][t] * backward[i][t];
    double denom = 0;

    for (int j = 0; j < numStates; j++) denom += forward[j][t] * backward[j][t];

    return divide(num, denom);
  }

  /** Prints all the parameters of an HiddenMarkovModelCSmoothed */
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
        System.out.print(
            "State(" + i + "," + j + ") = " + fmt.format(transitionMatrix[i][j]) + "  ");
      System.out.println();
    }

    System.out.println();
    System.out.println("Emission Matrix");
    emissionMatrix.forEach((k, v) -> System.out.println("Key: " + k + " Value: " + v));
  }

  /** divides two doubles. 0 / 0 = 0! */
  private double divide(double n, double d) {
    if (n == 0) return 0;
    else return n / d;
  }
}
