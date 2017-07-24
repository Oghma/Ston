/**
 * This class implements the trainer for the Ston Stemmer
 *
 * @author Matteo Lisotto (matteo.lisotto@gmail.com)
 */
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.stream.Stream;
import java.util.stream.Collectors;
import java.util.List;

public class Trainer {
  public static void main(String args[]) {
    String usage =
        "Usage: Trainer topology smoothed prefix suffix sigma_file test_file [path_tosave_model]\n topology: the topology of the hmm\n smoothed: smooth the probabilities (true or false)\n prefix: number of prefix states\n suffix: number of suffix states\n sigma_file: file that contain sigma\n test_file: file that contain the observations\n path_to_save_model: file to save the trained model";

    if (args.length >= 1 && (args[0].equals("-h") || args[0].equals("--help")))
      System.out.println(usage);
    else if (args.length >= 6) {
      String topology = args[0].toUpperCase();
      boolean smooth = args[1].equalsIgnoreCase("true");
      int numPrefixes = Integer.parseInt(args[2]);
      int numSuffixes = Integer.parseInt(args[3]);
      List<Character> sigma = null;
      List<String> test = null;
      HiddenMarkovModel hmm;

      try (Stream<Character> stream =
          Files.lines(Paths.get(args[4])).map(i -> (Character) i.charAt(0))) {
        sigma = stream.collect(Collectors.toList());
      } catch (IOException e) {
        e.printStackTrace();
      }

      try (Stream<String> stream = Files.lines(Paths.get(args[5]))) {
        test = stream.collect(Collectors.toList());
      } catch (IOException e) {
        e.printStackTrace();
      }

      hmm = StonModels.build(topology, smooth, numPrefixes, numSuffixes, sigma);
      System.out.println("Start training the hmm");
      hmm.train(test);
      System.out.println("Hmm trained");

      if (args.length == 7) HiddenMarkovModel.save(hmm, args[6]);
      else HiddenMarkovModel.save(hmm);
      System.out.println("Model saved");

    } else System.out.println(usage);
  }
}
