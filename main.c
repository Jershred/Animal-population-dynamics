/*Programme du projet 3
 ecrit par Jeremy Archier (p2019441)
  le 06/01/2021
  a=9
  b=4
  c=4
  d=1*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void coefficients(double t, double M, double *alphaN, double *alphaM, double *betaN, double *betaM)
{ /* Fonction qui calcule les coefficients alpha_i; beta_i */

    *betaN=-4.0/3.0; // Parametres de predation
    *betaM=1.0;

    if (t-floor(t)>=0.2 && t-floor(t)<0.5) // Condition d'introduction saisonniere de l'espece N
    {
        *alphaN=4.0/3.0;
    }
    else
    {
        *alphaN=2.0/3.0;
    }
    if (M>3) // Condition de chasse regulatrice de l'espece M
    {
        *alphaM=-1.5;
    }
    else
    {
        *alphaM=-1.0;
    }
}

void sauvegarde(double N, double M, double p, double t)
{ /* Fonction qui enregistre les valeurs de t, N(t), M(t), p(t) dans un fichier "population.dat" */

    FILE *fichier = NULL;  // On cree une variable qui pointe vers le fichier
    fichier = fopen("population.dat", "a");  // On va ecrire a la fin du fichier et le creer s'il n'existe pas encore
    fprintf(fichier, "t=%lf : N(t)=%lf  M(t)=%lf  p(t)=%lf \n",t,N,M,p);
    fclose(fichier);
}

void analyse_dynamique( double *N, double *M, int n)
{ /*fonction qui affiche a l'utilisateur le max de M et la moyenne de N*/

    double max=M[0];
    double somme=N[0];
    double moyenne;
    int i;

    for (i=1; i<n; i++)
    {
        if (M[i]>max)  // On determine le maximum de M
        {
            max=M[i];
        }
        somme=somme+N[i]; // On calcul la sommme de tous les valeurs de N
    }
    moyenne=somme/(double)n; // On obtient la moyenne pour n points
    printf("Le maximum de M est : %lf \n",max);
    printf("La moyenne de N est : %lf \n",moyenne);
}


void integre(double *alphaN, double *alphaM, double *betaN, double *betaM, int n)
{ /* Fonction qui integre le systeme d'equation (1),(2) */

    double *N=malloc(n*sizeof(double)); // Allocations memoire pour les differents tableaux
    double *M=malloc(n*sizeof(double));
    double *p=malloc(n*sizeof(double));
    double *t=malloc(n*sizeof(double));

    int i;
    double h=50.0/(double)(n-1); // h la longueur des sous-intervalles pour n points donc n-1 intervalles

    t[0]=0.0;  // Conditions initiales
    M[0]=2.0;
    N[0]=2.0;

    for (i=0; i<n; i++)
    {
        coefficients(t[i], M[i], alphaN, alphaM, betaN, betaM);  // On calcul les alpha i
        p[i]=M[i]/(M[i]+N[i]);  // Proportion de l'espece M
        sauvegarde(N[i], M[i], p[i], t[i]);  // On rentre les valeurs de t,N(t),M(t),p(t) dans le fichier "population.dat"

        N[i+1]=(1.0+*betaN*M[i]*h)/(1.0-*alphaN*h)*N[i]; // Il suffit d'utilise le theoreme des accroissements finis pour discretiser les termes de gauches du syteme de Lokta-Volterra
        M[i+1]=(1.0+*betaM*N[i]*h)/(1.0-*alphaM*h)*M[i]; // et ensuite en considerant les termes non-lineaires explicite et les termes lineaires implicite on obtient 2 relations de recurrence
        t[i+1]=t[i]+h;
    }

    analyse_dynamique(N, M, n); // On va afficher la dynamique de ces populations sur l'ensemble du cycle

    free(N); // On libere la memoire
    free(M);
    free(t);
    free(p);
}

int main()
{
    double alphaN, alphaM;  // On initialise les alpha i
    double betaN, betaM;  // On intitialise les beta i

    int n;  // Taille du tableau que l'on demande a l'utilisateur

    printf("Quelle taille de tableau desirez-vous pour les donnees t,N,M (il faut que la taille de tableau n>5000 pour que l'algorihtme converge) ? \n - ");
    // Le schéma d'Euler convergerait a partir de 100 points par seconde selon Loïc Reynier
    scanf("%d",&n);

    printf(" .\n .\n .\n Soyez patient.e le calcul peut prendre un peu longtemps en fonction du nombre de points choisis ;) \n .\n .\n .\n");

    integre(&alphaN, &alphaM, &betaN, &betaM, n); // On appelle la fonction integre qui appelera toutes les autres fonctions
    return 0;
}
