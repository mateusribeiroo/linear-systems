#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <stdbool.h>
//#include<stdint.h>

void freeMatriz(float **mat, int nlines)
{
    for(int i = 0; i < nlines; i++)
    {
        free(mat[i]);
    }
    free(mat);
}

void printMatriz(float **mat, int nlines, int mcols)
{
    for(int i = 0; i < nlines; i++)
    {
        for(int j = 0; j < mcols; j++)
        {
            printf("%.3f ", mat[i][j]);
            if(j == 2)
            {
                printf("| ");
            }
        }
        printf("\n");
    }
}

//TODO Fazer testes nessa funcao
bool sasselfeldCriterion(float **matH, int nlines, int mcols){
    float aux = 0, maior = 0;
    float beta[nlines];

    for(int i = 0; i < nlines; i++){
        beta[i] = 1;
    }

    for (int i = 0; i < nlines; i++)
    {
        for (int j = 0; j < mcols; j++)
        {
            beta[i] += matH[i][j]*beta[j];
        }
    }

    for(int i = 0; i < nlines; i++){
        if (i == 0) {
            beta[i] = maior;
        }
        if (beta[i] > maior) {
            maior = beta[i];
        }
    }

    if(maior > 1){
        return true;
    }

    return false;
}

bool errorVerification(float *vet1, float *vet2, int nlines, float erro){
    float maior = 0;
    float teste[nlines];

    for(int i = 0; i < nlines; i++)
    {
        teste[i] = fabs(vet2[i] - vet1[i])/fabs(vet2[i]);
    }


    for (int i = 0; i < nlines; i++)
    {
        if (i == 0) 
        {
            maior = teste[i];
        }

        if (teste[i] > maior) 
        {
            maior = teste[i];
        }
    }

    if(maior < erro)
    {
        printf("\nPrecisao minima atingida\n");
        return false;
    }

    return true;
}

void gaussSeidel(float **matH, float **mat, int nlines, int mcols, float erro, int max_it, float *pivos){
    int debug = 0, k = 0;
    float result[mcols], result2[mcols], vetAux[mcols];
    bool flag;

    printf("Para modo de Debug digite 1, para apenas ver os resultados sem os passos digite 0: ");
    scanf("%d", &debug);

    //testando criterio de sassenfeld
    if(sasselfeldCriterion(matH, nlines, mcols) == true){
        printf("Pelo critério de Sasselfeld o método de Gauss-Seidel não convergirá para esta matriz\n");
        return;
    }

    //Criando vetores de resultados
    for(int i = 0; i < mcols; i++)
    {
        if(i == mcols-1)
        {
            result[i] = 1;
            result2[i] = 1;
            vetAux[i] = 1;
            continue;
        }

        result[i] = 0;
        result2[i] = 0;
        vetAux[i] = 0;
    }
    

    //iniciar iteracoes---------------------------------------------------------
    float teste[mcols];
    printf("Solucoes de cada iteracao de Gauss-Seidel: \n");   
    do{
        for(int i = 0; i < nlines; i++)
        {
            for(int j = 0; j < mcols; j++)
            {
                vetAux[i] += matH[i][j] * result2[j];
            }

            result[i] = result2[i];
            result2[i] = vetAux[i];
            vetAux[i] = 0;
        }

        if(debug){
            printf("\n------------------------ Iteracao %d: --------------------------\n", k+1);
            for (int i = 0; i < nlines; i++)
            {
                printf("x(%d): %.3f\n", i+1, result2[i]);
            }
            printf("\n--------------------------------------------------------------\n");
        }

        k++;

        flag = errorVerification(result, result2, nlines, erro);
        
    }while(flag == true && k < max_it-1);

    printf("\n----------------------------------------------\n");
    printf("A solucao apromixada para este sistema é: ");
    for (int i = 0; i < nlines; i++)
    {
        printf("\nX(%d): %.5f ", i+1, result2[i]);
    }
    printf("\n----------------------------------------------\n");
}

bool linesCriterion(float **mat, int nlines, int mcols)
{
  return true;
}

void jacobiRichardson(float **matH, float **mat, int nlines, int mcols, float erro, int max_it, float *pivos)
{
    bool flag = true;
    if(linesCriterion(mat, nlines, mcols) != true)  
    {
        printf("A matriz NAO obedece ao criterio das linhas");
        return;
    }

    float result2[mcols], result[mcols], vetAux[mcols];
    
    //Criando vetores de resultados
    for(int i = 0; i < mcols; i++)
    {
        if(i == mcols-1)
        {
            result[i] = 1;
            result2[i] = 1;
            vetAux[i] = 1;
            continue;
        }

        result[i] = 0;
        result2[i] = 0;
        vetAux[i] = 0;
    }

    // iniciando iteracoes

    int k = 0;
    do {
        for (int i = 0; i < nlines; i++) 
        {
            for (int j = 0; j < mcols; j++) 
            {
                vetAux[i] += matH[i][j] * result[j];
            }

            result2[i] = vetAux[i];
            vetAux[i] = 0;
        }

        flag = errorVerification(result, result2, nlines, erro);

        for(int i = 0; i < nlines; i++)
        {
           result[i] = result2[i];
        }

        k++;
    }while (flag == true && k < max_it-1);

    printf("A solucao por Jacobi-Richardson: ");
    for(int i = 0; i < nlines; i++)
    {
        printf("\nX(%d): %.5f", i+1, result2[i]);
    }
    
    return;
}

void strictlyDominantDiagonal(float **mat, int nlines, int mcols, float *pivos)
{
    float soma[nlines];

    for(int i = 0; i < nlines; i++)
    {
        for(int j = 0; j < mcols; j++)
        {
            if(i != j && j != mcols-1)
            {
                soma[i] += fabs(mat[i][j]);
            }
        }

        if(soma[i] < pivos[i])
        {
            printf("\nA matriz NAO e estritamente diagonal dominante\n");
            return;
        }
    }
    
    printf("\nA matriz e SIM estritamente diagonal dominante\n");
}

void stantardLine(float **mat, int nlines, int mcols)
{
    int k = 0;
    float soma[nlines], maior = 0;

    for(int i = 0; i < nlines; i++)
    {
        for(int j = 0; j < mcols; j++)
        {
            if(i != j && j != mcols-1)
            {
                soma[i] += fabs(mat[i][j]);
            }
        }
    }

    for(int i = 0; i < nlines; i++)
    {
        if(i == 0)
        {
            maior = soma[i];
            k = i;
        }

        if(soma[i] > maior)
        {
            maior = soma[i];
            k = i;
        }
    }

    printf("\nA %dº e a norma-linha\n", k);
}


int main()
{
    float **mat, aux, maior, erro;
    int nlines, mcols, max_it, k = 0, choice;
    bool flag = true;

    while (flag == true)
    {
        printf("------------------------------------ Considere que voce esta criando uma matriz aumentada ------------------------------------\n");
        printf("Insira a quantidade de linhas da matriz: ");
        scanf("%d", &nlines);
        printf("Insira a quantidade de colunas da matriz: ");
        scanf("%d", &mcols);

        if(mcols - nlines != 1)
        {
            printf("\nA matriz precisa ser quadrada\n");
            continue;
        }

        printf("Insira a quantidade máxima de iteracoes que deseja calcular: ");
        scanf("%d", &max_it);

        if(max_it < 1)
        {
            printf("\nInsira 1 ou mais iteracoes\n");
            continue;
        }

        printf("Insira o erro maximo: ");
        scanf("%f", &erro);

        flag = false;
    }
    

    //alocando pointeiros do tipo float para cada posicao do vetor
    mat = malloc(nlines * sizeof (float*));

    //alocando um vetor de float para cada 
    for(int i = 0; i < nlines; i++)
    {
        mat[i] = malloc(mcols * sizeof(float));
    }

    float pivos[nlines];

    //iniciando a matriz com zeros para evitar erros de lixo de memória
    for(int i = 0; i < nlines; i++)
    {
        for(int j = 0; j < mcols; j++)
        {
            mat[i][j] = 0;
        }
    }

    printf("Veja sua matriz iniciada abaixo\n");
    printMatriz(mat, nlines, mcols);
    
    //recebe as entradas da matriz
    for(int i = 0; i < nlines; i++)
    {
        for(int j = 0; j < mcols; j++)
        {
            printf("\nInsira o valora para a posicao mat[%d][%d]: ", i+1, j+1);
            scanf("%f", &aux);
            mat[i][j] = aux;

            printMatriz(mat, nlines, mcols);

            if(i == j)
            {
                if(mat[i][j] == 0)
                {
                    printf("\nA diagonal principal nao pode ser nula\n");
                    freeMatriz(mat, nlines);
                    return 1;
                }
                pivos[i] = mat[i][j];
            }
        }
    }

    printf("\n------------------------------------------------------------------------\n");

    do{
        printf("\n------------------------------------ MENU ------------------------------------\n");
        printf("Digite 1 para usar o método de Gauss-Seidel\n");
        printf("Digite 2 para usar o método de Jacobi Richardson\n");
        printf("Digite 3 para ver a norma linha da sua matriz\n");
        printf("Digite 4 para descobrir se sua matriz é estritamente diagonal dominante\n");
        printf("Digite 5 para ver sua matriz H\n");
        printf("Digite 0 para sair do programa\n");
        scanf("%d", &choice);
        printf("------------------------------------------------------------------------\n");

        float **matH, result[mcols], result2[mcols], vetAux[mcols];
        matH = malloc(nlines * sizeof(float *));
        for(int i = 0; i < nlines; i++)
        {
            matH[i] = malloc(mcols * sizeof(float));
        }

        // calcula matriz H
        for (int i = 0; i < nlines; i++)
        {
            for (int j = 0; j < mcols; j++)
            {
                if(i == j){
                    matH[i][j] = 0;
                    continue;
                }

                if(j == mcols-1){
                    matH[i][j] = mat[i][j]/pivos[i];
                    continue;
                }

                if(j != mcols-1){
                    matH[i][j] = (-1 * mat[i][j])/pivos[i];
                    continue;
                }
            }
        }

        switch (choice)
        {
            case 1:
                gaussSeidel(matH, mat, nlines, mcols, erro, max_it, pivos);
                break;
            case 2:
                jacobiRichardson(matH, mat, nlines, mcols, erro, max_it, pivos);
                break;
            case 3:
                stantardLine(mat, nlines, mcols);
                break;
            case 4:
                strictlyDominantDiagonal(mat, nlines, mcols, pivos);
                break;
            case 5:
                printMatriz(matH, nlines, mcols);
                break;
            case 0:
                freeMatriz(mat, nlines);
                freeMatriz(matH, nlines);
                break;
            default:
                break;
        }

    }while(choice != 0);

    return 0;
}
