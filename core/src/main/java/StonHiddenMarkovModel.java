/**
 * This class implements the abstract class StonHiddenMarkovModel for the Ston Stemmer
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

public abstract class StonHiddenMarkovModel implements HiddenMarkovModel {
  /** Number of states */
  protected int numStates;

  /** Number of prefix states */
  protected int numPrefixes;

  /** Number of suffix states */
  protected int numSuffixes;

  /** Length of sigma */
  protected int sigmaSize;

  /** Sigma */
  protected List<Character> sigma;

  /** Initial state probabilities */
  protected double initialProbabilities[];

  /** Transition probabilities */
  protected double transitionMatrix[][];

  /** Emission probabilities */
  protected Map<Pair<Integer, Character>, Double> emissionMatrix;

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
        initialP[i] = divide(initialP[i], observations.size());
      }

      /* Re-estimation of transition probabilities */
      for (int i = 0; i < numStates; i++) {
        for (int j = 0; j < numStates; j++) {
          double num = 0;
          double denom = 0;
          for (int k = 0; k < numObservations; k++) {
            String observation = observations.get(k);
            for (int t = 0; t <= observation.length() - 1; t++) {
              num += epsilon(t, i, j, observation, forward.get(k), backward.get(k));
              denom += gamma(i, t, forward.get(k), backward.get(k));
            }
          }
          transitionM[i][j] = divide(num, denom);
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
          emissionM.put(new Pair<>(i, sigma.get(k)), divide(num, denom));
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
    int psi[][] = new int[numStates][lenObs];
    Integer path[] = new Integer[lenObs];
    double maxProb;
    int stemPosition;

    /* Calculate delta at time 0. At first finds the maxProb */
    for (int i = 0; i < numStates; i++) {
      delta[i][0] =
          initialProbabilities[i]
              * emissionMatrix.get(new Pair<Integer, Character>(i, observation.charAt(0)));
      psi[i][0] = 0;
    }

    /* Induction */
    for (int t = 1; t < lenObs; t++) {
      for (int j = 0; j < numStates; j++) {
        maxProb = 0;
        for (int k = 0; k < numStates; k++) {
          double prob = delta[k][t - 1] * transitionMatrix[k][j];
          if (prob > maxProb) {
            maxProb = prob;
            psi[j][t] = k;
          }
        }
        delta[j][t] =
            maxProb * emissionMatrix.get(new Pair<Integer, Character>(j, observation.charAt(t)));
      }
    }

    /* Termination */
    path[lenObs - 1] = maxFinalState(delta, lenObs);

    /* Calculates the path backtracking */
    for (int t = lenObs - 2; t >= 0; t--) {
      path[t] = psi[path[t + 1]][t + 1];
    }

    stemPosition = Arrays.asList(path).indexOf(numPrefixes);
    return stemPosition != -1 ? observation.substring(0, stemPosition) : observation;
  }

  /**
   * Calculation of Forward-Variables f(i,t) for state i at time t for output sequence O with the
   * current HiddenMarkovModelC parameters
   *
   * @param observations the output sequence O
   * @return an array f(i,t) over states and times, containing the Forward-variables.
   */
  private List<Double[][]> forwardProc(List<String> observations) {
    List<Double[][]> forward = new ArrayList<Double[][]>();

    for (String observation : observations) {
      int lenObs = observation.length();
      Double forwardI[][] = new Double[numStates][lenObs];

      /* Initialization (time 0) */
      for (int i = 0; i < numStates; i++)
        forwardI[i][0] =
            initialProbabilities[i] * emissionMatrix.get(new Pair<>(i, observation.charAt(0)));

      /* Induction */
      for (int t = 0; t <= lenObs - 2; t++) {
        for (int j = 0; j < numStates; j++) {
          forwardI[j][t + 1] = 0.0;
          for (int i = 0; i < numStates; i++)
            forwardI[j][t + 1] += (forwardI[i][t] * transitionMatrix[i][j]);
          forwardI[j][t + 1] *= emissionMatrix.get(new Pair<>(j, observation.charAt(t + 1)));
        }
      }
      forward.add(forwardI);
    }

    return forward;
  }

  /**
   * Calculation of Backward-Variables b(i,t) for state i at time t for output sequence O with the
   * current HiddenMarkovModelC parameters
   *
   * @param observations the output sequence O
   * @return an array b(i,t) over states and times, containing the Backward-Variables.
   */
  private List<Double[][]> backwardProc(List<String> observations) {
    List<Double[][]> backward = new ArrayList<Double[][]>();

    for (String observation : observations) {
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
                    * emissionMatrix.get(new Pair<>(j, observation.charAt(t + 1))));
        }
      }
      backward.add(backwardI);
    }

    return backward;
  }

  /**
   * Calculation of probability P(X_t = s_i, X_t+1 = s_j | O, m).
   *
   * @param t time t
   * @param i the number of state s_i
   * @param j the number of state s_j
   * @param observable an output sequence o
   * @param forward the Forward-Variables
   * @param backward the Backward-Variables
   * @return epsilon of the state s_i
   */
  private double epsilon(
      int t, int i, int j, String observable, Double[][] forward, Double[][] backward) {
    double num;
    if (t == observable.length() - 1) num = forward[i][t] * transitionMatrix[i][j];
    else
      num =
          forward[i][t]
              * transitionMatrix[i][j]
              * emissionMatrix.get(new Pair<>(j, observable.charAt(t + 1)))
              * backward[j][t + 1];

    double denom = 0;

    for (int k = 0; k < numStates; k++) denom += (forward[k][t] * backward[k][t]);

    return divide(num, denom);
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

  abstract int maxFinalState(double delta[][], int lenObs);
}
