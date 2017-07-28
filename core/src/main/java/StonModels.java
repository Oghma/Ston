import java.util.List;
import java.util.Map;
import javafx.util.Pair;

public class StonModels {
  public static HiddenMarkovModel build(
      String topology, boolean smoothed, int numPrefixes, int numSuffixes, List<Character> sigma) {
    HiddenMarkovModel hmm = null;
    switch (topology) {
      case "A":
        hmm =
            smoothed == true
                ? new HiddenMarkovModelASmoothed(numPrefixes, numSuffixes, sigma)
                : new HiddenMarkovModelA(numPrefixes, numSuffixes, sigma);
        break;
      case "B":
        hmm =
            smoothed == true
                ? new HiddenMarkovModelBSmoothed(numPrefixes, numSuffixes, sigma)
                : new HiddenMarkovModelB(numPrefixes, numSuffixes, sigma);

        break;
      case "C":
        hmm =
            smoothed == true
                ? new HiddenMarkovModelCSmoothed(numPrefixes, numSuffixes, sigma)
                : new HiddenMarkovModelC(numPrefixes, numSuffixes, sigma);
        break;
    }
    return hmm;
  }
}
