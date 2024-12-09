# Equation-de-la-chaleur-1D_2D
Ce projet consiste en la simulation de l'equation de la chaleur
en 1D et en 2D, avec une resolution numerique avec les differences finies et des conditions de Dirichlet.

On utilisera la bibliothèque `numpy` pour la résolution et `vpython` pour la visualisation.

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

#### 1D

#### 2D

# Méthode implocite 

#### 1D

#### 2D


