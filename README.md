# Interpolation Numérique : Lagrange, Newton et Hermite

Ce projet étudie et compare trois méthodes d'interpolation polynomiale (**Lagrange**, **Newton** et **Hermite**) appliquées à deux fonctions aux comportements très différents sur l'intervalle $[-4, 4]$ :
1. La fonction exponentielle 
2. La fonction de Runge 

##  Structure du Projet

- **`interpolation.py`** : Contient l'implémentation des algorithmes de calcul des coefficients :
  - `computePL` : Formulation de Lagrange.
  - `compute_pN` : Différences divisées de Newton.
  - `compute_pH` : Interpolation d'Hermite (prenant en compte les valeurs et les dérivées).
- **`main.py`** : Script principal qui génère le maillage, évalue les polynômes et trace la grille comparative de 6 graphiques ($2 \times 3$).

---

## Analyse des Résultats et Comparaison des Graphiques

### 1. Équivalence entre Lagrange et Newton
- **Observation :** Les tracés pour Lagrange et Newton sont **strictement identiques** sur tous les graphiques.

### 2. Comportement sur la fonction exponentielle 
- **Observation :** À mesure que le nombre de nœuds $N$ augmente ($N = 3, 5, 9, 13$),
 l'erreur d'interpolation diminue rapidement sur tout l'intervalle $[-4, 4]$.

### 3. Le Phénomène de Runge sur la fct rationnelle
- **Observation :** Au centre de l'intervalle ($x \approx 0$),
 l'interpolation devient très précise quand $N$ augmente. En revanche, **près des bords ($x \approx -4$ et $x \approx 4$), 
 les polynômes développent de fortes oscillations** qui s'amplifient considérablement pour $N = 9$ et $N = 13$.

### 4. Apport de l'interpolation d'Hermite
- **Observation :** Hermite offre une meilleure approximation du sommet de la courbe (en $x=0$)
 dès les plus petites valeurs de $N$.
