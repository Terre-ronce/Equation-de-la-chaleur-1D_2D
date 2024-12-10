from vpython import *
import numpy as np
import bisect

diffusivite_thermique_materiaux = { #en mm2/s
    "Acier": 22.8,
    "Aluminium": 98.8,
    "Argent": 172,
    "Calcium": 123,
    "Cuivre": 117,
    "Diamant": 306,
    "Eau liquide": 0.144,
    "Eau solide": 0.541,
    "Fer": 20.4,
    "Gold": 129,
    "Nickel": 22.3,
    "Or": 129,
    "Plomb": 23.3,
    "Tin": 37,
    "Tungstène": 76.1,
    "Verre": 0.543,
    "Zinc": 40.2
}

def temperature_a_couleur(temp, temp_min, temp_max):
    """
    Renvoie une couleur en fonction de la température normée:
    - 0.0 à 0.2 : Bleu foncé → Bleu
    - 0.2 à 0.4 : Bleu → Vert turquoise
    - 0.4 à 0.6 : Vert turquoise → Jaune clair
    - 0.6 à 0.8 : Jaune clair → Orange
    - 0.8 à 1.0 : Orange → Rouge → Rouge foncé
    
    Arguments:
    - temp : float - Température actuelle
    - temp_min : float - Température minimale de la série
    - temp_max : float - Température maximale de la série
    
    Retourne:
    - vector : Couleur VPython correspondante.
    """
    # on normalise la température entre [0, 1]
    norm_temp = (temp - temp_min) / (temp_max - temp_min)
    norm_temp = max(0, min(1, norm_temp))  # précaution

    # on initialise les couleurs
    r, g, b = 0, 0, 0

    # on transcrit les températures normés en couleurs alant du bleu foncé au rouge foncé en
    # passant par le jaune, le vert et le orange ainsi que des nuances plus clair de ces couleurs
    if norm_temp <= 0.2:
        r = 0
        g = 0
        b = 0.5 + 2.5 * norm_temp
    elif norm_temp <= 0.4:
        r = 0
        g = (norm_temp - 0.2) / 0.2
        b = 1
    elif norm_temp <= 0.6:
        r = (norm_temp - 0.4) / 0.2
        g = 1
        b = 1 - r
    elif norm_temp <= 0.8:
        r = 1
        g = 1 - (norm_temp - 0.6) / 0.4
        b = 0
    else:
        r = 1 - (norm_temp - 0.8) / 0.4
        g = 0
        b = 0

    return vector(r, g, b)

def matrice_diagonale(taille, valeur):
    """
    Génère une matrice diagonale de taille n × n.
    
    Arguments :
    taille : int - Taille de la matrice (n × n).
    valeur : float - Valeur sur la diagonale principale.
    
    Retourne :
    numpy.ndarray - Matrice diagonale de taille n × n.
    """
    # on nitialise la matrice n × n avec des zéros
    M = np.zeros((taille, taille))
    
    # on rempli la diagonale principale avec la valeur
    np.fill_diagonal(M, valeur)

    return M

def matrice_tridiagonale(n, c1, c2):
    """
    Génère une matrice tridiagonale de taille n × n.
    
    Arguments :
    n : int - Taille de la matrice (n × n).
    c1 : float - Valeur sur la diagonale principale.
    c2 : float - Valeur sur les diagonales immédiatement au-dessus et au-dessous.
    
    Retourne :
    numpy.ndarray - Matrice tridiagonale de taille n × n.
    """
    # on nitialise la matrice n × n avec des zéros
    M = np.zeros((n, n))
    
    # on rempli la diagonale principale avec la valeur
    np.fill_diagonal(M, c1)
    
    # on rempli les diagonales supérieure et inférieure avec la valeur
    np.fill_diagonal(M[:-1, 1:], c2)  # supérieure
    np.fill_diagonal(M[1:, :-1], c2)  # inférieure
    
    return M

def resoudre_equation_matricielle(A, B):
    """
    Résout l'équation matricielle AX = B.
    
    Arguments :
    A : numpy.ndarray - Matrice carrée de coefficients (taille n × n).
    B : numpy.ndarray - Matrice ou vecteur des termes constants.
    
    Retourne :
    numpy.ndarray - Solution X de l'équation AX = B.
    
    Lève une erreur si A est singulière (non inversible).
    """
    try:
        # résout l'équation AX = B
        X = np.linalg.solve(A, B)
        return X
    except np.linalg.LinAlgError as e:
        raise ValueError("La matrice A est singulière ou mal conditionnée.") from e

def trouver_encadrement(valeur, liste):
    """
    Trouve les indices des deux éléments de la liste triée qui encadrent la valeur donnée.
    
    Arguments:
    - valeur : float - La valeur à encadrer.
    - liste : list[float] - Une liste triée de flottants.
    
    Retourne:
    - tuple (int, int) : Les indices des deux éléments qui encadrent la valeur.
                         Si la valeur est en dehors de la liste, retourne (None, None) ou des indices ajustés.
    """
    # on trouve l'indice d'insertion par dichotomie
    index = bisect.bisect_left(liste, valeur)
    
    # si l'élément se placerait à la place du premier élément
    if index == 0:
        return None, 0
    # si l'élément se placerait à la place du dernier élément
    elif index == len(liste):
        return len(liste) - 1, None
    
    # on retourne les indices encadrants
    return index - 1, index


def moyenne_ponderee(index1, index2, y, x, valeur):
    """
    Calcule la moyenne pondérée des valeurs y[index1] et y[index2]
    en utilisant une valeur intermédiaire pour définir les poids basés sur les distances.

    Arguments:
    - index1 : int - Premier index (à gauche).
    - index2 : int - Second index (à droite).
    - y : list[float] - Liste des valeurs associées (y).
    - x : list[float] - Liste des positions ou distances (x).
    - valeur : float - Valeur pour laquelle on calcule la moyenne pondérée.

    Retourne:
    - float : La moyenne pondérée.
    """
    # vérification de cohérences
    if index1 < 0 or index2 >= len(x) or index1 >= index2:
        raise ValueError("Indices invalides ou hors des limites.")
    if not (x[index1] <= valeur <= x[index2]):
        raise ValueError("La valeur n'est pas encadrée par les indices donnés.")
    
    # on calcule les poids relatifs
    distance_total = x[index2] - x[index1]
    poids1 = (x[index2] - valeur) / distance_total
    poids2 = (valeur - x[index1]) / distance_total

    # on retourne la moyenne pondérée
    return poids1 * y[index1] + poids2 * y[index2]