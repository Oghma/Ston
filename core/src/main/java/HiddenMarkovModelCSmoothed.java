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

public class HiddenMarkovModelCSmoothed extends StonHiddenMarkovModelSmoothed {
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
