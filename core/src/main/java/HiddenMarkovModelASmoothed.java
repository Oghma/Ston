/**
 * This class implements the topology A of Hidden Markov Model for the Ston Stemmer
 *
 * @author Matteo Lisotto (matteo.lisotto@gmail.com)
 */
import java.util.List;
import java.util.Map;
import javafx.util.Pair;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.Arrays;
import java.text.DecimalFormat;

public class HiddenMarkovModelASmoothed extends StonHiddenMarkovModelSmoothed {
  /** Contain the final states of the chain */
  private List<Integer> finalStates;

  /**
   * Initializes an HiddenMarkovModelASmoothed.
   *
   * @param numPrefixes number of prefix states
   * @param numSuffixes number of suffix states
   * @param sigma the alphabet
   */
  public HiddenMarkovModelASmoothed(int numPrefixes, int numSuffixes, List<Character> sigma) {
    this.numPrefixes = 2;
    for (int i = 1; i <= numPrefixes - 2; i++) this.numPrefixes += i;

    this.numSuffixes = 0;
    for (int i = 1; i <= numSuffixes; i++) this.numSuffixes += i;

    this.numStates = this.numPrefixes + this.numSuffixes;

    this.initialProbabilities = new double[numStates];
    this.transitionMatrix = new double[numStates][numStates];
    this.emissionMatrix = new LinkedHashMap<>();
    this.sigma = new ArrayList<>(sigma);
    this.smooth = new ArrayList<Double[]>();
    this.finalStates = new ArrayList<>();

    int state = 0;
    int lineStates = 1;
    double initialProbability = 1.0 / (numPrefixes - 2);

    for (int i = 0; i < numPrefixes - 2; i++) {
      /* Initializes initalProbabilities parameter */
      initialProbabilities[state] = initialProbability;

      /* Initializes the prefixes in the transitionMatrix parameter  */
      double level[][] = createInitialLevel(lineStates);
      copyRegion(level, state);
      state += lineStates;

      transitionMatrix[state - 1][this.numPrefixes - 2] = 1.0;

      lineStates++;
    }
    transitionMatrix[0][this.numPrefixes - 2] = 0.5;
    transitionMatrix[this.numPrefixes - 2][this.numPrefixes - 1] = 1.0;

    /* Initializes the suffixes in the transitionMatrix parameter and finalState */
    finalStates.add(this.numPrefixes - 1);
    state = this.numPrefixes;
    lineStates = 1;
    for (int i = 0; i < numSuffixes; i++) {
      double level[][] = createFinalLevel(lineStates);
      copyRegion(level, state);
      transitionMatrix[this.numPrefixes - 1][state] = 1.0;

      state += lineStates;
      lineStates++;
      finalStates.add(state - 1);
    }

    /* Initializes the emission probabilities */
    this.sigmaSize = sigma.size();
    double emissionProbability = 1.0 / sigmaSize;

    for (int i = 0; i < numStates; i++) {
      for (Character c : sigma)
        emissionMatrix.put(new Pair<Integer, Character>(i, c), emissionProbability);
    }
  }

  private double[][] createInitialLevel(int numStates) {
    double states[][] = new double[numStates][numStates];
    states[0][0] = 0.5;
    for (int i = 1; i < numStates - 1; i++) {
      states[i][i + 1] = 1.0;
    }

    if (numStates > 1) states[0][1] = 0.5;

    return states;
  }

  private double[][] createFinalLevel(int numStates) {
    double states[][] = new double[numStates][numStates];
    for (int i = 0; i < numStates - 1; i++) {
      states[i][i + 1] = 1.0;
    }
    return states;
  }

  private void copyRegion(double[][] level, int state) {
    int lenLevel = level[0].length;

    for (int i = 0; i < lenLevel; i++) {
      for (int j = 0; j < lenLevel; j++) transitionMatrix[state + i][state + j] = level[i][j];
    }
  }

  int maxFinalState(double delta[][], int lenObs) {
    double maxProb = 0;
    int maxState = 0;

    for (int state : finalStates)
      if (delta[state][lenObs - 1] > maxProb) {
        maxProb = delta[state][lenObs - 1];
        maxState = state;
      }
    return maxState;
  }
}
