/**
 * This class implements the stemmer for the Ston Stemmer
 *
 * @author Matteo Lisotto (matteo.lisotto@gmail.com)
 */
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.stream.Stream;
import java.util.stream.Collectors;
import java.util.List;

public class Stemmer {
  public static void main(String args[]) {
    String usage =
        "Usage: Stemmer stem_file [model_file] [stemmed_file]\n stem_file: file that contain the words to stem\n model_file: file that contain the hmm\n stemmed_file: where to save stemmed words";

    if (args.length >= 1) {
      if (args[0].equals("-h") || args[0].equals("--help")) System.out.println(usage);
      else {
        List<String> stem = null;
        HiddenMarkovModel hmm;
        List<String> stemmed;

        try (Stream<String> stream = Files.lines(Paths.get(args[0]))) {
          stem = stream.collect(Collectors.toList());
        } catch (IOException e) {
          e.printStackTrace();
        }

        hmm = args.length >= 2 ? HiddenMarkovModel.load(args[1]) : HiddenMarkovModel.load();
        System.out.println("Hmm loaded");

        System.out.println("Starting decoding");
        stemmed = stem.stream().map(w -> hmm.decode(w)).collect(Collectors.toList());
        System.out.println("Words stemmed");

        try {
          if (args.length == 3) Files.write(Paths.get(args[2]), stemmed);
          else Files.write(Paths.get("stemmedWords.txt"), stemmed);
        } catch (IOException e) {
          e.printStackTrace();
        }
        System.out.println("Stemmed words saved");
      }
    } else System.out.println(usage);
  }
}
