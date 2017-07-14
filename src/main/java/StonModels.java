import java.util.Map;
import javafx.util.Pair;

public class StonModels {
    public static HiddenMarkovModel build(String topology, int numPrefixes, int numSuffixes) {
	HiddenMarkovModel hmm = null;
	switch(topology) {
	case "A":
	    break;
	case "B":
	    break;
	case "C":
	    hmm = new HiddenMarkovModelC(numPrefixes, numSuffixes);
	    break;
	}
	return hmm;
    }

    public static HiddenMarkovModel build(String topology, int numPrefixes, int numSuffixes, Map<Pair<Integer, Character>, Double> emissions) {
	HiddenMarkovModel hmm = null;
	switch(topology) {
	case "A":
	    break;
	case "B":
	    break;
	case "C":
	    hmm = new HiddenMarkovModelC(numPrefixes, numSuffixes, emissions);
	    break;
	}
	return hmm;
    }
}
