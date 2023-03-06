import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;
import java.lang.Integer;

public class GenerateInput {
    public static void generateAndWriteAdjacencyMatrixToFile(int n, double p, String filename) {
        int[][] matrix = new int[n][n];
        Random random = new Random();

        // Preenche a matriz com 0s e 1s aleatórios
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    matrix[i][j] = 0;
                } else {
                    if (random.nextDouble() < p) {
                        Random w = new Random();
                        matrix[i][j] = w.nextInt(15) + 1;
                        matrix[j][i] = w.nextInt(15) + 1;
                    } else {
                        matrix[i][j] = 99999;
                        matrix[j][i] = 99999;
                    }
                }
            }
        }

        // Escreve a matriz em um arquivo
        try {
            BufferedWriter writer = new BufferedWriter(new FileWriter(filename));
            writer.write(n + " ");
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    writer.write(matrix[i][j] + " ");
                    if (j == n - 1)
                        writer.write("\n");
                }
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

   public static void main(String[] args){
        generateAndWriteAdjacencyMatrixToFile(Integer.valueOf(args[0]), 0.85, "../matrix.txt");
    }
}