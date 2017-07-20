import java.util.List;
import java.io.FileOutputStream;
import java.io.ObjectOutputStream;
import java.io.ObjectOutput;
import java.io.IOException;

public interface HiddenMarkovModel implements Serializable {

  public void train(List<String> data);

  public void train(List<String> data, int steps);

  public void print();

  public String decode(String observation);

  public static void save(HiddenMarkovModel hmm) {
    try (FileOutputStream outputFile = new FileOutputStream("model.data");
        ObjectOutputStream saveModel = new ObjectOutputStream(outputFile)) {
      saveModel.writeObject(hmm);
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

  public static void save(HiddenMarkovModel hmm, String path) {
    try (FileOutputStream outputFile = new FileOutputStream("model.data");
        ObjectOutputStream saveModel = new ObjectOutputStream(outputFile)) {
      saveModel.writeObject(hmm);
    } catch (IOException e) {
      e.printStackTrace();
    }
  }
}
