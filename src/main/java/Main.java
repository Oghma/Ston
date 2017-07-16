import java.util.*;

import javafx.util.Pair;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class Main {
    public static void main (String args[]) {
	List<String> test = new ArrayList<>();
	List<String> testSet = new ArrayList<>();
	String fileName = "";
    int[] array = new Random().ints(50, 0, 386609).toArray();
    ArrayList<Character> sigma
                = new ArrayList<>(
                "qwertyuioplkjhgfdsazxcvbnm1234567890".chars()
                        .mapToObj(c -> (char) c)
                        .collect(Collectors.toList()));


    try (Stream<String> stream = Files.lines(Paths.get(fileName))) {
                test = stream.collect(Collectors.toList());
	} catch (IOException e) {
	    e.printStackTrace();
	}

	for (int i = 0; i < array.length; i++) {
	    int j = array[i];
	    String s = test.get(j);
	    testSet.add(s);
    }



	
	HiddenMarkovModel hmm = StonModels.build("C", 7,3, sigma);
	System.out.println(testSet);
	System.out.println("Training");
	hmm.train(testSet, 1);
	hmm.print();
    }

    static Map<Pair<Integer,Character>, Double>  buildEmissions() {
	Map<Pair<Integer,Character>, Double> em = new LinkedHashMap<>();
	String alp = "qwertyuioplkjhgfdsazxcvbnm1234567890";
	double prob = 1.0 / alp.length();

	
	for(int i = 0; i < 10; i++) {
	    for(int j = 0; j < alp.length(); j++) {
		em.put(new Pair<>(i, alp.charAt(j)), prob);
	    }
	}

	return em;
    }
}
