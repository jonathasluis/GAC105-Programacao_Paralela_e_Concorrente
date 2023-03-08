# GAC105-Programacao_Paralela_e_Concorrente
Repositório com arquivos GAC105 - Programação Paralela e Concorrente - UFLA 2022/2 

## Trabalho de implementação MPI

### Modo de Usar: 
- Passo 1:
  - Compilar o "GenerateInput.java" no diretório "Input".
  - Executar o programa Java passando como parâmetro a quantidade de Vértices.
  - Isso vai gerar o arquivo "matrix.txt" que representa a Matriz de Adjacência de um grafo.

- Passo 2:
  - Compilar o arquivo "Main.cpp".
  - Executar o programa gerado.
    - Parametros:
      - "-c" compila o "Bellmanford.cpp";
      - "-c -r" compila e roda "Bellmanford.cpp";
      - Se nenhum parameatro for passado, apenas roda "Bellmanford.cpp"
      - OBS: ao rodar "Bellmanford.cpp" será executado 30 vezes, com a mesma quantidade de execuções para cada quantidade de processos (1, 2 e 4)
  - A execução desses Programas gerarão arquivos texto na pasta "output"

- Passo 3:
    - O arquivo "timeMedias.txt" contém as informações das medidas, basta apenas importa-lo na planilha.
