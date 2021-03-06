/**
 * This class implements the abstract class HiddenMarkovModelSmoothed for the Ston Stemmer
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

public abstract class StonHiddenMarkovModelSmoothed implements HiddenMarkovModel {
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

  /** Smoothing variable */
  protected List<Double[]> smooth;

  /** Serial version UID generated with serialver */
  private static final long serialVersionUID = -383347686438850352L;
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
          emissionMatrix.put(new Pair<Integer, Character>(i, sigma.get(k)), divide(num, denom));
        }
      }
      initialProbabilities = initialP;
      transitionMatrix = transitionM;
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

    /* Calculate delta at time 0. */
    for (int i = 0; i < numStates; i++) {
      delta[i][0] =
          Math.log(initialProbabilities[i])
              + Math.log(
                  emissionMatrix.get(new Pair<Integer, Character>(i, observation.charAt(0))));
      psi[i][0] = 0;
    }

    /* Induction */
    for (int t = 1; t < lenObs; t++) {
      for (int j = 0; j < numStates; j++) {
        maxProb = Double.NEGATIVE_INFINITY;
        for (int k = 0; k < numStates; k++) {
          double prob = delta[k][t - 1] + Math.log(transitionMatrix[k][j]);
          if (prob > maxProb) {
            maxProb = prob;
            psi[j][t] = k;
          }
        }
        delta[j][t] =
            maxProb
                + Math.log(
                    emissionMatrix.get(new Pair<Integer, Character>(j, observation.charAt(t))));
      }
    }

    /* Termination */
    maxProb = Double.NEGATIVE_INFINITY;
    for (int i = 0; i < numStates; i++)
      if (delta[i][lenObs - 1] > maxProb) {
        maxProb = delta[i][lenObs - 1];
        path[lenObs - 1] = i;
      }

    /* Calculates the path backtracking */
    for (int t = lenObs - 2; t >= 0; t--) {
      path[t] = psi[path[t + 1]][t + 1];
    }

    stemPosition = Arrays.asList(path).indexOf(numPrefixes - 1);
    return stemPosition != -1 ? observation.substring(0, stemPosition) : observation;
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

  abstract int maxFinalState(double delta[][], int lenObs);
}
