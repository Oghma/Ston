import java.util.List;
import java.util.Map;
import javafx.util.Pair;

public class StonModels {
    public static HiddenMarkovModel build(String topology, int numPrefixes, int numSuffixes, List<Character> sigma) {
	HiddenMarkovModel hmm = null;
	switch(topology) {
	case "A":
	    break;
	case "B":
	    break;
	case "C":
	    hmm = new HiddenMarkovModelC(numPrefixes, numSuffixes, sigma);
	    break;
	}
	return hmm;
    }

    public static HiddenMarkovModel build(String topology, int numPrefixes, int numSuffixes, List<Character> sigma, Map<Pair<Integer, Character>, Double> emissions) {
	HiddenMarkovModel hmm = null;
	switch(topology) {
	case "A":
	    break;
	case "B":
	    break;
	case "C":
	    hmm = new HiddenMarkovModelC(numPrefixes, numSuffixes, sigma, emissions);
	    break;
	}
	return hmm;
    }
}
