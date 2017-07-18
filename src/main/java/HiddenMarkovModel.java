import java.util.List;

public interface HiddenMarkovModel {

  public void train(List<String> data);

  public void train(List<String> data, int steps);

  public void print();

  public String decode(String observation);
}
