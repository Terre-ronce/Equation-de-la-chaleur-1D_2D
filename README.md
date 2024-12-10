# Equation-de-la-chaleur-1D_2D
Ce projet consiste en la simulation de l'équation de la chaleur
en 1D (barre à section rectangulaire) et en 2D (plaque), avec une résolution numérique avec les differences finies et des conditions de Dirichlet.

On utilisera la bibliothèque `numpy` pour la résolution et `vpython` pour la visualisation.

### Restrictions
Les simulations se font dans les limites suivantes :

- Conditions aux limites de Dirichlet (comprendre que les conditions aux limites sont maintenues le long de la simulation, il n'y a pas de variation) ;
- Pour les plaques on prendra des plaques carrés avec un même pas spatiale suivant les deux directions ;

### Exemple

#### 1D

Pour le problème à une dimension les paramètres sont :

- Dimension : $1D$ ;
- Matériau : $Aluminium$ ;
- Diffusivité thermique : $98.8 mm^2.s^{-1}$ ;
- Longueur : $75 mm$ ;
- Nombre de points d'espaces : $40$ ;
- Nombre de points de temps : $1071$ ;
- Durée : $4 s$.
- Température initiale dans la barre : $20°C$ ;
- Condition limite à gauche : $90°C$ ;
- Condition limite à droite : $90°C$.

On parvient à la simulation suivante (vidéo accélérée) :

<p align="center">
  <img src="https://github.com/user-attachments/assets/eebaca7a-4c8e-4125-afe5-0e6f511b2a5c" alt="animation" />
</p>

#### 2D

De même en deux dimensions, en conservant les paramètre précédents et en ajoutant :

- Hauteur : $75 mm$ ;
- Nombre de points d'espaces en $x$: $40$ ;
- Nombre de points d'espaces en $y$: $40$ ;
- Durée : $4 s$.
- Condition limite en haut : $90°C$ ;
- Condition limite en bas : $90°C$.

On parvient à la simulation suivante (vidéo accélérée) :

<p align="center">
  <img src="https://github.com/user-attachments/assets/03ae5f2b-6d19-4a5f-ae5e-951a91cb01b8" alt="animation" />
</p>

## Equation de la chaleur

L'équation générale de la chaleur sans source de chaleur est donnée par :

$$
\frac{\partial u}{\partial t} = D \nabla^2 u
$$

où :
- $u = u(x, y, z, t)$ est la température en fonction de l'espace et du temps,
- $D$ est la diffusivité thermique (constante positive),
- $\nabla^2$ est l'opérateur laplacien.

### Cas en 1D

En une dimension spatiale $(x)$, l'équation devient :

$$
\frac{\partial u}{\partial t} = \alpha \frac{\partial^2 u}{\partial x^2}
$$

### Cas en 2D

En deux dimensions spatiales $(x, y)$, l'équation devient :

$$
\frac{\partial u}{\partial t} = \alpha \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right)
$$

## Résolution numérique
### Méthode explicite 

On peut résoudre numériquement l'équation de la chaleur de manière explicite enécrivant les différences finies avec les instants $n$ et $n+1$ pour les dérivés temporelles et l'instant n, pour les dérivés spatiales [[1]](https://www.f-legrand.fr/scidoc/docmml/numerique/diffusion/methode/diffusion.html). En un position donnée, on peut alors calculé la température à l'instant $n+1$ connaissant la température à l'instant $n$. 


#### 1D

En réécrivant l'équation 1D de la manière suivante [[2]](https://cahier-de-prepa.fr/psi-encpb/download?id=687) :

$$
\frac{u_j^{n+1} - u_j^{n}}{\Delta t} = D \frac{u_{j-1}^{n} - 2u_{j}^{n} + u_{j+1}^{n}}{\Delta x^2}
$$

On a donc la formule explicite :

$$
u_j^{n+1} = u_j^{n} + r (u_{j-1}^{n} - 2u_{j}^{n} + u_{j+1}^{n})
$$

Avec $r = D\frac{\Delta t}{\Delta tx^2}$.
#### 2D

De même, en 2D on a [[3]](https://www.f-legrand.fr/scidoc/docmml/numerique/diffusion/diffusion2d/diffusion2d.html) :

$$
u_{i,j}^{n+1} = u_{i,j}^{n} + r (u_{i-1,j}^{n} - 2u_{i,j}^{n} + u_{i+1,j}^{n}) + r (u_{i,j-1}^{n} - 2u_{j}^{n} + u_{i,j+1}^{n})
$$

Avec **$r = D\frac{\Delta t}{\Delta x^2} = D\frac{\Delta t}{\Delta y^2}$**.

#### Observation
Ce schéma est précis est au premier ordre et stable à la condition $r\leq0.5$. Condition qui peut imposer un pas de temps et/ou un pas d'espace trop élevé.

# Méthode implicite 
Le schéma implicite avec le modèle de Crank-Nicholson, est précis au second ordre et stable sans condition (mais plus gourmand en calcul) [1].

#### 1D
En 1D et dans les mêmes notations que précedemment on a la formule implicite suivante [2] :

$$
u_{j+1}^n = u_j^n + r^2 \big( u_j^{n+1} - 2u_j^n + u_j^{n-1} + u_{j+1}^{n+1} - 2u_{j+1}^n + u_{j+1}^{n-1} \big)
$$

#### 2D
De même, on a [[4]](https://www.ijser.org/researchpaper/A-Case-study-on-simulation-of-heat-equation-by-Crank-Nicolson-Method-in-Accordance-with-digital-image-processing.pdf) :

$$
u_{j+1,k}^n = u_{j,k}^n + \frac{r}{2} \big[ 
    \big( u_{j+1,k}^{n+1} - 2u_{j,k}^{n+1} + u_{j-1,k}^{n+1} + u_{j,k+1}^{n+1} - 2u_{j,k}^{n+1} + u_{j,k-1}^{n+1} \big) + 
    \big( u_{j+1,k}^n - 2u_{j,k}^n + u_{j-1,k}^n + u_{j,k+1}^n - 2u_{j,k}^n + u_{j,k-1}^n \big)
\big]
$$

#### Observation
En pratique on écrira le scéma implicite sous forme matricielle. Les détails des formules se trouvent dans la littérature [2][4].

## Simulation

Après la résolution on crée des corps discrédités (succession de petits parallélogrammes), on assigne des couleurs à chaque valeur de température et enfin pour chaque instant de simulation, on met à jour les couleurs pour pouvoir observer l'évolution.

Le programme choisi lui même en fonction de la valeur de $r$, quel schéma adopter pour la résolution (mais il est possible de forcer un schéma).

## Limites

Je réalise pour mes exemples des calculs sur de petites durée et dimensions n'ayant pas beaucoup de puissance avec ma machine, mais il est possible pour ceux ayant des machines plus performantes de faire plus conséquent. Toutefois, le fait que le programme soit en python reste un frein.

## Références

[1] : Frédéric Legrand. Équation de diffusion à une dimension.

[2] : La PSI de l'ENCPB. Méthode des différences finies : propagation thermique dans une poutre.

[3] : Frédéric Legrand. Équation de diffusion à deux dimensions.

[4] : Irfan Raju et al. A Case study on simulation of heat equation by Crank-Nicolson Method in Accordance with digital image processing
