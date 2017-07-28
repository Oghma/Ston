/**
 * This class implements the topology B of Hidden Markov Model for the Ston Stemmer
 *
 * @author Matteo Lisotto (matteo.lisotto@gmail.com)
 */
import java.util.Arrays;
import java.util.LinkedHashMap;
import javafx.util.Pair;
import java.util.List;
import java.util.ArrayList;

public class HiddenMarkovModelB extends StonHiddenMarkovModel {
  /** Contain the final states of the chain */
  private List<Integer> finalStates;

  /** Serial version UID generated with serialver */
  private static final long serialVersionUID = 1470784236577942433L;

  /**
   * Initializes an HiddenMarkovModelB.
   *
   * @param numPrefixes number of prefix states
   * @param numSuffixes number of suffix states
   * @param sigma the alphabet
   */
  public HiddenMarkovModelB(int numPrefixes, int numSuffixes, List<Character> sigma) {
    this.numPrefixes = numPrefixes;
    this.numSuffixes = numSuffixes;

    this.numSuffixes = 0;
    for (int i = 1; i <= numSuffixes; i++) this.numSuffixes += i;

    this.numStates = this.numPrefixes + this.numSuffixes;

    this.initialProbabilities = new double[numStates];
    this.transitionMatrix = new double[numStates][numStates];
    this.emissionMatrix = new LinkedHashMap<>();
    this.sigma = new ArrayList<>(sigma);
    this.finalStates = new ArrayList<>();

    /* Initialize the parameter initialProbability */
    double initialProbability = 1.0 / (numPrefixes - 2);
    Arrays.fill(this.initialProbabilities, 0, numPrefixes - 2, initialProbability);

    for (int i = 0; i < numPrefixes - 2; i++) {
      this.transitionMatrix[i][i] = 0.5;
      this.transitionMatrix[i][i + 1] = 0.5;
    }

    for (int i = numPrefixes - 2; i < numPrefixes - 1; i++) {
      this.transitionMatrix[i][i + 1] = 1;
    }

    /* Initializes the suffixes in the transitionMatrix parameter and finalState */
    finalStates.add(this.numPrefixes - 1);
    int state = this.numPrefixes;
    int lineStates = 1;
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

    for (int i = numPrefixes - 1; i < numStates; i++)
      if (delta[i][lenObs - 1] > maxProb) {
        maxProb = delta[i][lenObs - 1];
        maxState = i;
      }
    return maxState;
  }
}
