import java.util.List;
import java.io.FileOutputStream;
import java.io.ObjectOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.Serializable;
import java.io.FileInputStream;

public interface HiddenMarkovModel extends Serializable {

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

  public static HiddenMarkovModel load() {
    HiddenMarkovModel hmm = null;
    try (FileInputStream inputFile = new FileInputStream("model.data");
        ObjectInputStream loadModel = new ObjectInputStream(inputFile)) {
      hmm = (HiddenMarkovModel) loadModel.readObject();
    } catch (IOException i) {
      i.printStackTrace();
    } catch (ClassNotFoundException c) {
      System.out.println("HiddenMarkovModel class not found");
      c.printStackTrace();
    }
    return hmm;
  }

  public static HiddenMarkovModel load(String path) {
    HiddenMarkovModel hmm = null;
    try (FileInputStream inputFile = new FileInputStream(path);
        ObjectInputStream loadModel = new ObjectInputStream(inputFile)) {
      hmm = (HiddenMarkovModel) loadModel.readObject();
    } catch (IOException i) {
      i.printStackTrace();
    } catch (ClassNotFoundException c) {
      System.out.println("HiddenMarkovModel class not found");
      c.printStackTrace();
    }
    return hmm;
  }
}
