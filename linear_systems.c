#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <stdbool.h>
//#include<inttypes>

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
int sasselfeldCriterion(float **matH, int nlines, int mcols){
    float aux = 0, maior = 0;
    float beta[nlines];

    for(int i = 0; i < nlines; i++){
        beta[i] = 1;
    }

    for (int i = 0; i < nlines; i++)
    {
        for (int j = 0; j < mcols; j++)
        {
            if(i == 0){
                beta[i] += matH[i][j];
            }else{
                aux += matH[i][j]*beta[i];
            }
        }
        beta[i] = aux;
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
        return 1;
    }

    return 0;
}

bool errorVerification(float *vet1, float *vet2, int size, float erro){
    float maior = 0;
    float teste[size];

    for(int i = 0; i < size-1; i++)
    {
        printf("\n%.3f - %.3f/ %.3f\n", vet2[i], vet1[i], vet2[i]);
        teste[i] = fabs(vet2[i] - vet1[i])/fabs(vet2[i]);
    }

    for (int i = 0; i < size; i++)
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
        return false;
    }

    return true;
}

void gaussSeidel(float **matH2, float **matH, float **mat, int nlines, int mcols, float erro, int max_it, float *pivos){
    int debug = 0, k = 0;
    float result[mcols], result2[mcols], vetAux[mcols];
    bool flag;

    printf("Para modo de Debug digite 1, para apenas ver os resultados sem os passos digite 0: ");
    scanf("%d", &debug);

    //testando criterio de sassenfeld
    if(sasselfeldCriterion(matH, nlines, mcols)){
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

    // matriz H recebe os valores da matriz base
    for (int i = 0; i < nlines; i++)
    {
        for (int j = 0; j < mcols; j++)
        {
            matH[i][j] = mat[i][j];
        }
    }

    // calcula matriz H
    for (int i = 0; i < nlines; i++)
    {
        for (int j = 0; j < mcols; j++)
        {
            if(i == j){
                matH2[i][j] = 0;
                continue;
            }

            if(j == mcols-1){
                matH2[i][j] = matH[i][j]/pivos[i];
                continue;
            }

            if(j != mcols-1){
                matH2[i][j] = (-1 * matH[i][j])/pivos[i];
                continue;
            }
        }
    }

    if(debug == 1){
        printf("\nSua matriz H: ------------------------------------------------------\n");
        printMatriz(matH2, nlines, mcols);
        printf("------------------------------------------------------\n");
    }

    //iniciar iteracoes---------------------------------------------------------
    float teste[mcols];
    printf("Solucoes de cada iteracao de Gauss-Seidel: \n");   
    do{
        for(int i = 0; i < nlines; i++)
        {
            for(int j = 0; j < mcols; j++)
            {
                result2[i] += matH2[i][j] * vetAux[j];
            }

            result[i] = vetAux[i];
            vetAux[i] = result2[i];
            result2[i] = 0;
        }

        if(debug == 1){
            printf("\n------------------------Iteracao %d: --------------------------\n", k+1);
            for (int i = 0; i < mcols-1; i++)
            {
                printf("x(%d): %.3f\n", i+1, vetAux[i]);
            }
            printf("\n--------------------------------------------------------------\n");
        }

        k++;

        flag = errorVerification(result, vetAux, mcols, erro);

        for (int i = 0; i < mcols-1; i++)
            result[i] = result2[i];
        
    }while(flag == true && k < max_it-1);

    printf("\n----------------------------------------------\n");
    printf("A solucao apromixada para este sistema é: ");
    for (int i = 0; i < mcols-1; i++)
    {
        printf("\nX(%d): %.3f ", i+1, vetAux[i]);
    }
    printf("\n----------------------------------------------\n");
}

//TODO
void jacobiRichardson(float **matH2, float **matH, float **mat, int nlines, int mcols, float erro, int max_it, float *pivos)
{

}

void strictlyDominantDiagonal(float **mat, float *pivos)
{

}

void stantardLine(float **mat)
{

}


int main()
{
    float **mat, aux, maior, erro;
    int nlines, mcols, max_it, k = 0;
    int choice;
    bool flag;

    printf("------------------------------------ Considere que voce esta criando uma matriz aumentada ------------------------------------\n");
    printf("Insira a quantidade de linhas da matriz: ");
    scanf("%d", &nlines);
    printf("Insira a quantidade de colunas da matriz: ");
    scanf("%d", &mcols);
    printf("Insira a quantidade máxima de iteracoes que deseja calcular: ");
    scanf("%d", &max_it);
    printf("Insira o erro maximo: ");
    scanf("%f", &erro);

    //alocando pointeiros do tipo float para cada posicao do vetor
    mat = malloc(nlines * sizeof (float*));

    //alocando um vetor de float para cada 
    for(int i = 0; i < nlines; i++)
    {
        mat[i] = malloc(mcols * sizeof(float));
    }

    float pivos[nlines - 1];

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
        printf("Digite 0 para sair do programa\n");
        scanf("%d", &choice);
        printf("------------------------------------------------------------------------\n");

        float **matH, **matH2, result[mcols], result2[mcols], vetAux[mcols];
        matH = malloc(nlines * sizeof(float *));
        for(int i = 0; i < nlines; i++)
        {
            matH[i] = malloc(mcols * sizeof(float));
        }

        matH2 = malloc(nlines * sizeof(float *));
        for(int i = 0; i < nlines; i++)
        {
            matH2[i] = malloc(mcols * sizeof(float));
        }

        switch (choice)
        {
            case 1:
                gaussSeidel(matH2, matH, mat, nlines, mcols, erro, max_it, pivos);
                break;
            case 2:
                jacobiRichardson(matH2, matH, mat, nlines, mcols, erro, max_it, pivos);
                break;
            case 3:
                stantardLine(mat);
                break;
            case 4:
                strictlyDominantDiagonal(mat, pivos);
                break;
            case 0:
                freeMatriz(mat, nlines);
                freeMatriz(matH, nlines);
                freeMatriz(matH2, nlines);
                break;
            default:
                break;
        }

    }while(choice != 0);

    return 0;
}