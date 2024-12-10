from vpython import *
import numpy as np
import utils

class EquationChaleur:
    """
    Classe pour résoudre l'équation de la chaleur en 1D ou 2D et simuler la diffusion.
    Pour le cas 2D, on se placera dans le cas particulier d'une structure carrée (Nx = Ny et dx = dy).

    Attributs:
        dimension (str): La dimension de l'équation ('1D' ou '2D').
        materiaux (str): Le matériau de la barre.
        L (float): La longueur (et hauteur en 2D) de la barre.
        T (float): La température initiale en °C.
        Tg (float): La température initiale sur la face gauche, en °C.
                    Peut être une liste dans le cas 2D de deux éléments pour une température en triangle.
        Td (float): La température initiale sur la face droite, en °C.
                    Peut être une liste dans le cas 2D de deux éléments pour une température en triangle.
        Tb (float): La température initiale sur la face basse, en °C.
                    Peut être une liste dans le cas 2D de deux éléments pour une température en triangle.
        Th (float): La température initiale sur la face haute, en °C.
                    Peut être une liste dans le cas 2D de deux éléments pour une température en triangle.
        duree (float): La durée de la simulation.
        Nt (int): Le nombre de pas de temps.
        Nx (int): Le nombre de pas d'espace.
        dt (float): Le pas de temps.
        dx (float): Le pas d'espace (x ou y).
        r (float): Un coefficient.
        corps (list): Les éléments graphiques de la simulation.
        centres (list): Les centres des éléments graphiques.
        Tfs (list): Les températures pour chaque point de la barre à chaque instant.
        vitesse (int): La vitesse de simulation (nombre de mise à jour par seconde).

    Méthodes:
        EulerExplicite(): Résout l'équation de la chaleur en utilisant la méthode d'Euler explicite.
        EulerImplicite(): Résout l'équation de la chaleur en utilisant la méthode d'Euler implicite avec Crank-Nicholson.
        ObtenirCouleur(): Calcule les couleurs correspondant aux températures pour chaque point de la simulation.
        CreerElement(): Crée la représentation graphique 3D de la barre ou de la plaque.
        Simuler(): Lance la simulation de la diffusion thermique, permet aussi de suivre l'évolution de la température en
                    un point donné de l'élément.
    """

    def __init__(self,dimension: str, materiaux: str, L: float = 50, T: float = 20,
                duree: float = 4, Nt: int = 1000, Nx: int = 40):
        """
        Initialise la classe EquationChaleur.

        Arguments:
            dimension (str): La dimension de l'équation ('1D' ou '2D').
            materiaux (str): Le matériau de la barre.
            L (float): La longueur de la barre.
            T (float): La température initiale.
            duree (float): La durée de la simulation.
            Nt (int): Le nombre de pas de temps.
            Nx (int): Le nombre de pas d'espace.
        """

        self.dimension = dimension
        self.materiaux =materiaux
        # diffusivité thermique
        self.D = utils.diffusivite_thermique_materiaux[self.materiaux]

        self.Nt = Nt
        self.Nx =Nx

        self.T = T
        # nombre différents de conditions aux limites (Dirichlet) différents en 1D et 2D.
        if self.dimension == '1D':
            self.Tg = input('Température à gauche (nombre décimal (en °C)) : ')
            self.Td = input('Température à droite (nombre décimal (en °C)) : ')
            # vitesse pour une simulation en 1D
            self.vitesse = 1000
        elif self.dimension == '2D':
            self.Tb = eval(input('Température en bas (nombre décimal ou liste de deux éléments (en °C)) : '))
            self.Th = eval(input('Température en haut (nombre décimal ou liste de deux éléments (en °C)) : '))
            self.Tg = eval(input('Température à gauche (nombre décimal ou liste de deux éléments (en °C)) : '))
            self.Td = eval(input('Température à droite (nombre décimal ou liste de deux éléments (en °C)) : '))
            # vitesse pour une simulation en 2D
            self.vitesse = 10000
        else:
            raise ValueError("`dimension` doit être `1D` ou `2D`.")
        
        self.L = L
        self.e = 8

        self.duree = duree

        self.dt = self.duree / self.Nt

        self.dx = self.L / (self.Nx - 1)
        self.r= self.D * self.dt / self.dx ** 2

        self.corps = []
        self.centres = [ - self.L/2 + i * self.dx + self.dx/2 for i in range(self.Nx)]

        # températures finales on utilise la résolution stable la moins couteuse.
        self.Tfs = self.EulerExplicite() if self.r < 0.5 else self.EulerImplicite()

    def EulerExplicite(self):
        """
        Résolution de l'équation de la chaleur en utilisant la méthode d'Euler explicite.
        """

        # en 1D
        if self.dimension == '1D':
            # température initiale dans la barre
            finale = [np.zeros(self.Nx) + self.T]
            finale[0][0] = self.Tg
            finale[0][-1] = self.Td

            # à chaque instant
            for _ in range(self.Nt - 1):
                # températures dans la barre à l'instant d'avant
                T_avant = finale[-1]
                # initialisation pour températures à cet instant
                T_mtn = np.zeros(self.Nx)
                # on les calcule avec la formule
                for i in range(1, self.Nx - 1):
                    T_mtn[i] = T_avant[i] + self.r * (T_avant[i - 1] - 2 * T_avant[i] + T_avant[i + 1])
                # on réinitialise les températures aux bords (car CL Dirichlet)
                T_mtn[0] = self.Tg
                T_mtn[-1] = self.Td
                # on enregistre les températures pour cet instant et on passe à l'instant suivant
                finale.append(T_mtn)
        # en 2D
        elif self.dimension == '2D':
            # températures initiales dans la plaque
            finale = [np.zeros((self.Nx, self.Nx)) + self.T]
            # conditions aux limites
            finale[0][0, :] = np.linspace(self.Td[0], self.Td[1], self.Nx) if isinstance(self.Td, list) else self.Td
            finale[0][-1, :] = np.linspace(self.Tg[0], self.Tg[1], self.Nx) if isinstance(self.Tg, list) else self.Tg
            finale[0][:, 0] = np.linspace(self.Th[0], self.Th[1], self.Nx) if isinstance(self.Th, list) else self.Th
            finale[0][:, 1] = np.linspace(self.Tb[0], self.Tb[1], self.Nx) if isinstance(self.Tb, list) else self.Tb

            # à chaque instant
            for _ in range(self.Nt - 1):
                # températures dans la plaque à l'instant d'avant
                T_avant = finale[-1]
                # initialisation pour températures à cet instant
                T_mtn = np.zeros_like(T_avant)
                # on les calcule avec la formule
                for i in range(1, self.Nx - 1):
                    for j in range(1, self.Nx - 1):
                        T_mtn[i, j] = (
                            T_avant[i, j]
                            + self.r * (T_avant[i - 1, j] - 2 * T_avant[i, j] + T_avant[i + 1, j])
                            + self.r * (T_avant[i, j - 1] - 2 * T_avant[i, j] + T_avant[i, j + 1])
                        )

                # on réinitialise les températures aux bords (car CL Dirichlet)
                T_mtn[0, :] = np.linspace(self.Td[0], self.Td[1], self.Nx) if isinstance(self.Td, list) else self.Td
                T_mtn[-1, :] = np.linspace(self.Tg[0], self.Tg[1], self.Nx) if isinstance(self.Tg, list) else self.Tg
                T_mtn[:, 0] = np.linspace(self.Th[0], self.Th[1], self.Nx) if isinstance(self.Th, list) else self.Th
                T_mtn[:, -1] = np.linspace(self.Tb[0], self.Tb[1], self.Nx) if isinstance(self.Tb, list) else self.Tb
                # on enregistre les températures pour cet instant et on passe à l'instant suivant
                finale.append(T_mtn)

        return finale

    def EulerImplicite(self):
        """
        Résout l'équation de la chaleur en utilisant la méthode d'Euler implicite.
        On utilise Crank-Nicholson pour plus de robustesse et de stabilité.
        """

        # en 1D
        if self.dimension == '1D':
            # température initiale dans la barre
            finale = [np.zeros(self.Nx) + self.T]
            # conditions aux limites
            finale[0][0] = self.Tg
            finale[0][-1] = self.Td

            # matrices de l'équation
            A = utils.matrice_tridiagonale(self.Nx, 1 + self.r, -self.r/ 2)
            M = utils.matrice_tridiagonale(self.Nx, 1 - self.r, self.r/ 2)
            A[0][0], A[self.Nx-1][self.Nx-1] = 1, 1
            A[0][1], A[self.Nx-1][self.Nx-2] = 0, 0
            M[0][0], M[self.Nx-1][self.Nx-1] = 1, 1
            M[0][1], M[self.Nx-1][self.Nx-2] = 0, 0

            # à chaque instant
            for _ in range(self.Nt - 1):
                # on résout l'équation matricielle pour avoir les températures
                B = np.dot(M, finale[-1])
                # on enregistre 
                finale.append(utils.resoudre_equation_matricielle(A, B))
                # on réinitialise les conditions aux limites (car CL Dirichlet) et on passe à l'instant suivant
                finale[-1][0] = self.Tg
                finale[-1][-1] = self.Td
        
        # en 2D
        elif self.dimension == '2D':
            # températures initiales dans la plaque
            finale = [np.zeros((self.Nx, self.Nx)) + self.T]
            # conditions aux limites
            finale[0][0, :] = np.linspace(self.Td[0], self.Td[1], self.Nx) if isinstance(self.Td, list) else self.Td
            finale[0][-1, :] = np.linspace(self.Tg[0], self.Tg[1], self.Nx) if isinstance(self.Tg, list) else self.Tg
            finale[0][:, 0] = np.linspace(self.Th[0], self.Th[1], self.Nx) if isinstance(self.Th, list) else self.Th
            finale[0][:, -1] = np.linspace(self.Tb[0], self.Tb[1], self.Nx) if isinstance(self.Tb, list) else self.Tb

            # matrices de l'équation
            A = utils.matrice_tridiagonale(self.Nx, 2 + 2 * self.r, -self.r)
            M = utils.matrice_tridiagonale(self.Nx, 2 - 2 * self.r, self.r)

            # à chaque instant
            for _ in range(self.Nt - 1):
                # températures en x
                T_col = np.zeros((self.Nx, self.Nx))
                # on résout
                for j in range(self.Nx):
                    B = np.dot(M, finale[-1][:, j])
                    T_col[:, j] = utils.resoudre_equation_matricielle(A, B)
                # températures en y
                T_ligne = np.zeros((self.Nx, self.Nx))
                # on résous
                for i in range(self.Nx):
                    B = np.dot(M, T_col[i, :])
                    T_ligne[i, :] = utils.resoudre_equation_matricielle(A, B)

                # on réinitialise les températures aux bords (car CL Dirichlet)
                T_ligne[0, :] = np.linspace(self.Td[0], self.Td[1], self.Nx) if isinstance(self.Td, list) else self.Td
                T_ligne[-1, :] = np.linspace(self.Tg[0], self.Tg[1], self.Nx) if isinstance(self.Tg, list) else self.Tg
                T_ligne[:, 0] = np.linspace(self.Th[0], self.Th[1], self.Nx) if isinstance(self.Th, list) else self.Th
                T_ligne[:, -1] = np.linspace(self.Tb[0], self.Tb[1], self.Nx) if isinstance(self.Tb, list) else self.Tb

                # on enregistre et on passe à l'instant suivant
                finale.append(T_ligne)


        return finale

    @property
    def ObtenirCouleur(self):
        """
        Calcule les couleurs correspondant aux températures pour chaque point de la simulation
        et à chaque instant, en fonction de la dimension (1D ou 2D).
        """

        # en 1D
        if self.dimension == '1D':
            # on recupère les températures extrêmes initiales
            Tmin = min(min(self.Tfs[i]) for i in range(len(self.Tfs)))
            Tmax = max(max(self.Tfs[i]) for i in range(len(self.Tfs)))

            # on initialise
            couleurs = []
            # à chaque instant
            for i in range(len(self.Tfs)):
                # on calcule la couleur pour chaque point en fonction de sa température
                col = [
                    utils.temperature_a_couleur(float(self.Tfs[i][j]), Tmin, Tmax)
                    for j in range(len(self.Tfs[i]))
                ]
                # on enregistre
                couleurs.append(col)
        # en 2D
        elif self.dimension == '2D':
            # on recupère les températures extrêmes initiales
            Tmin = min(min(row) for t in self.Tfs for row in t)
            Tmax = max(max(row) for t in self.Tfs for row in t)

            # on initialise
            couleurs = []
            # à chaque instant
            for i in range(len(self.Tfs)):
                # on initialise
                temp_couleurs = []
                # pour chaque plaque
                for x in range(len(self.Tfs[i])):
                # on calcule la couleur pour chaque point en fonction de sa température
                    row_couleurs = [
                        utils.temperature_a_couleur(float(self.Tfs[i][x][y]), Tmin, Tmax)
                        for y in range(len(self.Tfs[i][x]))
                    ]
                    # on enregistre la colonne
                    temp_couleurs.append(row_couleurs)
                # on enregistre la plaque
                couleurs.append(temp_couleurs)

        return couleurs
    
    def CreerElement(self):
        """
        Crée la représentation graphique de la barre.
        """
        # en 1D
        if self.dimension == '1D':
            # pour chaque subdivision
            for i in range(self.Nx):
                    # on calcule la position du centre
                    position = - self.L/2 + i * self.dx + self.dx/2
                    # on crée le petit élément
                    portion = box(pos=vector(position, 0, 0), size=vector(self.dx, self.e, self.e), color=color.blue, emissive=True)
                    # on l'enregistre
                    self.corps.append(portion)

        # en 2D
        elif self.dimension == '2D':
            # pour chaque subdivision x
            for i in range(self.Nx):
                # on initialise
                rang = []
                # pour chaque subdivision y
                for j in range(self.Nx):
                    # on calcule la position du centre
                    position_x = -self.L / 2 + i * self.dx + self.dx / 2
                    position_y = -self.L / 2 + j * self.dx + self.dx / 2
                    # on crée le petit élément le long de x
                    bloc = box(pos=vector(position_x, position_y, 0),
                                size=vector(self.dx, self.dx, self.e),
                                color=color.blue, emissive=True)
                    # on l'enregistre
                    rang.append(bloc)
                # on enregistre la portion
                self.corps.append(rang)
    
    def Simuler(self, suivre_point=False, rapport=None):
        """
        Lance la simulation de la diffusion et, si demandé, suit un ou plusieurs points
        pour afficher leur évolution de température en temps réel.

        Arguments :
            suivre_point (bool) : Si True, suit un point pendant la simulation.
            rapport (float ou list[float]) :
                - Pour une barre (1D), un seul rapport (float) indique la position relative sur la longueur L.
                - Pour une plaque (2D), une liste [rapport_x, rapport_y] indique la position relative dans la plaque.
        """

        # on crée de l'élément (barre ou plaque)
        self.CreerElement()

        # on récupère les couleurs
        couleurs = self.ObtenirCouleur

        pos_x, pos_y = (- self.L / 7, - 3 * self.L / 5) if self.dimension == '2D' else (- 5 * self.e / 3, - self.e)

        # on initialise le temps de simulation restant
        temps_restant_label = label(pos=vector(pos_x, pos_y, 0),
                                    text="Temps restant pour la diffusion : {:.2f} s".format(self.duree),
                                    box=False, height=12, color=color.white)

        # suivi de point si demandé
        if suivre_point:
            temperatures = []

            if self.dimension == '1D':
                if rapport is None:
                    raise ValueError("Veuillez fournir un rapport pour le suivi du point en 1D.")
                # on calcul la position et indices encadrants
                point = self.L * rapport - self.L / 2
                index_1, index_2 = utils.trouver_encadrement(point, self.centres)
                titre = f"Évolution de la température au point d'abcisse {self.L * np.array(rapport)}."
            elif self.dimension == '2D':
                if rapport is None or not isinstance(rapport, list) or len(rapport) != 2:
                    raise ValueError("Veuillez fournir un rapport [rapport_x, rapport_y] pour le suivi du point en 2D.")
                # on calcul les positions et indices encadrants
                point_x = self.L * rapport[0] - self.L / 2
                point_y = self.L * rapport[1] - self.L / 2
                index_x1, index_x2 = utils.trouver_encadrement(point_x, self.centres)
                index_y1, index_y2 = utils.trouver_encadrement(point_y, self.centres)
                titre = f"Évolution de la température au point de coordonnée ({self.L * rapport[0]}, {self.L * rapport[1]})."

            # on crée le graphique
            graph(
                title=titre,
                xtitle="Temps (s)",
                ytitle="Température (°C)",
                width=600,
                height=400,
                background=color.white,
            )
            temperature_curve = gcurve(color=color.red)
        
        # simulation de la diffusion
        for t in range(self.Nt):
            # on met à jour le temps restant
            temps_ecoule = t * self.dt
            temps_restant = self.duree - temps_ecoule
            temps_restant_label.text = "Temps restant pour la diffusion : {:.2f} s".format(temps_restant)

            if self.dimension == '1D':
                for j in range(self.Nx):
                    rate(self.vitesse)
                    self.corps[j].color = couleurs[t][j]

                # on met à jour le suivi du point en 1D
                if suivre_point:
                    temp = utils.moyenne_ponderee(index_1, index_2, self.Tfs[t], self.centres, point)
                    temperatures.append(temp)
                    temperature_curve.plot(temps_ecoule, temp)

            elif self.dimension == '2D':
                for i in range(self.Nx):
                    for j in range(self.Nx):
                        rate(self.vitesse)
                        self.corps[i][j].color = couleurs[t][i][j]

                # on met à jour le suivi du point en 2D
                if suivre_point:
                    temp_x1 = utils.moyenne_ponderee(index_y1, index_y2, self.Tfs[t][index_x1], self.centres, point_y)
                    temp_x2 = utils.moyenne_ponderee(index_y1, index_y2, self.Tfs[t][index_x2], self.centres, point_y)
                    tempx = (temp_x1 + temp_x2)/2
                    temp_y1 = utils.moyenne_ponderee(index_x1, index_x2, self.Tfs[t][index_y1], self.centres, point_x)
                    temp_y2 = utils.moyenne_ponderee(index_x1, index_x2, self.Tfs[t][index_y2], self.centres, point_x)
                    tempy = (temp_y1 + temp_y2)/2
                    temp = (tempx + tempy)/2
                    temperatures.append(temp)
                    temperature_curve.plot(temps_ecoule, temp)

        # on force le temps à 0s à la fin de la simulation
        temps_restant_label.text = "Temps restant pour la diffusion : {:.2f} s".format(0.00)