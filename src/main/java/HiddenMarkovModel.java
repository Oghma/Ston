public interface HiddenMarkovModel {

    public <T> void train(T data);
    public <T> void train(T data, int steps);
    public void print();
    public String decode(String observation);
}
