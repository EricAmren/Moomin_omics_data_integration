# Introduction 

(présenter le contexte scientifique dans lequel s'inscrit le travail réalisé ainsi que la problématique abordée pendant le stage )

* La cellule, unité de base du vivant. Sans elle, vie impossible car pas d'équilibre : interface avec milieu extérieur, plus ou moins perméable selon la matière (et ses propriétés physico-chimiques...) et les conditions environnementales (pression, température etc).
* Maintien de l'intégrité de la cellule grâce à une régulation adaptée des échanges avec l'extérieur et entre ses différentes composantes
* Ces échanges se font par le biais de réaction chimique
* Réaction chimique = transformation de la matière
* Espèces consommées = réactifs
* Espèces formées  = produits
* réseau métabolique = ensemble des réactions chimiques d'une cellule
* Divisible en sous-ensemble fonctionnel : voie métabolique
* Les petites molécules élémentaires échangées = métabolites
* Perturbation de l'environnement (stress) induit la nécessité d'une réponse de la cellule pour rétablir son équilibre chimique
* Réponse par le biais d'un remaniement de son réseau métabolique
* => C'est ce mécanisme que l'on souhaiterait comprendre davantage : comment simplifier ce problème pour plus facilement l'appréhender ?
* Réseau métabolique modélisable sous la forme d'un graphe en joignant les réactions chimiques entre elles (via leurs réactifs et produits) 
* Graphe = ensemble de sommets et d'arêtes. ( métabolites et réactions chimiques de manière interchangeables selon les besoins de la modélisation)
* Réseau métabolique paramétrable via la régulation de la vitesse des réactions chimiques
* Biologiquement, plusieurs modalités peuvent influencer les vitesses de réaction (température, cofacteurs biochimiques)
* Mais ce sont majoritairement les enzymes qui déterminent la vitesse des réactions chimiques et par extension l'activité des voies métaboliques d'une cellule.
* Enzyme = protéine ayant la propriété de réduire l'énergie d'activation d'une réaction chimique et donc d'accélérer cette réaction : on dit qu'elle catalyse cette réaction.
* La production d'une enzyme peut être accrue (ou réduite) selon les conditions, et notamment selon l'abondance d'une molécule : c'est ce que l'on appelle l'induction (ou inhibition) enzymatique.
* Cette régulation se fait directement au niveau du taux d'expression du gène codant pour cette enzyme et c'est un des mécanismes principales de régulation du réseau métabolique.
* Pour rappel, un gène, constitué d'ADN, est transcrit en ARN : on dit qu'il est exprimé.
* Cette quantité d'ARN transcrite par gène est variable entre différentes conditions et ce que l'on cherche à mesurer lors d'une analyse d'expression différentielle.
* L'ARN est ensuite traduit en protéine.
* Certaines protéines peuvent alors jouer le rôle d'enzyme pour catalyser des réactions.
* Si l'on connaît les concentrations de ces enzymes avant et après une perturbation, on pourra avoir une idée des modifications opérées au sein du réseau métabolique.
* Il est donc possible de caractériser la réponse d'une cellule à un stress à plusieurs niveaux, d'amont en aval : (Donner avantages/inconvénients pour chaque) (peut-être un schéma ?)
  * Analyse du différentiel d'expression génique (mesure de l'expression des gènes codant pour des enzymes)
  * Analyse de données protéomiques (mesure des concentrations des enzymes)
  * Observer directement les flux (FBA) (mesure directe des échanges)
  * Analyse de données métabolomiques (mesure des différences de concentration en métabolite provoquées par ces changements)
* L'hypothèse mise en tension par mon laboratoire est qu'en tenant compte de résultats provenant de plusieurs méthodes, les conclusions n'en seront que plus robustes.
* Citer Kotoura, Totoro et Moomin
* Recentrer sur Moomin : pour l'instant, uniquement possible sur des données transcriptomiques provenant d'analyse RNASeq
* Il existe un autre type d'analyse de différentiel d'expression utilisant des puces à ADN. Cependant la sortie est différente et incompatible avec Moomin.
* Décrire différences entre ces deux types de données (comptages bruts VS ratios de luminosité, valeurs absolues VS relatives)
* => Problématique : comment intégrer des données de microarrays au logiciel moomin ? Est-il pertinent d'utiliser ce type de données ? Les prédictions obtenues par la suite sont-elles cohérentes avec celles obtenues avec des données RNASeq ?

# Matériels et méthodes

## Méthodes :

### Moomin :

* Entrée : 
  * matrice contenant les noms de gènes et leur fold-change
  * Modèle métabolique
* Sortie : génère des hypothèses de changements métaboliques 
* Algorithme : Problème d'optimisation linéaire 

### EBSeq

* méthode bayésienne permettant d'obtenir :
* Valeur de log fold-change pour chaque gène entre condition
* Valeurs de PPDE pour chaque gène
* EBSeq prend en entrée des données de comptage de transcrits (=données RNASeq)

## Jeux de données :

* Premier jeu de microarrays du papier 
* résultats de Limma
* Jeu de données RNAseq du papier :
  * Paired-end
  * nombre de réplicats  et de conditions
* résultats de FeatureCounts
* Deuxième jeu de microarrays pour la preuve de concept
* Modèles métaboliques (Core + IJO1366)



## Matériels :

* PC portable (2.40GHz * 4, 8Go de RAM)
* Pedago-NGS

# Résultats

## Comparaison des résultats de DEG obtenus par les différentes méthodes

* Intersection des résultats et proportions
* Tests statistiques :
  * Linear modelling (R²)
  * Spearman ($\rho​$)
* Graphes des logFC
* Volcano plots logFC/PPDE ou logFC/Pval

## Comparaison des reconstructions des modèles obtenus par Moomin selon les méthodes

* Intersection des résultats
* Tests statistiques :
  * Same same
* Petit exemple de cartes Escher

# Discussion

* Intérêt d'utiliser des méthodes bayésiennes (notamment dans les cas où l'on peut affirmer qu'il n'y a pas de différentiel d'expression) 

# Conclusion et perspectives : 

## Conclusion :

* => Au cours de ce stage, nous avons testé la faisabilité d'utiliser des données provenant de microarrays pour caractériser les remaniements opérés sur le réseau métabolique grâce au logiciel Moomin.
* Nous avons pour cela comparé les différentiels d'expression obtenus avec cette méthode avec ceux d'autres méthodes : EBSeq, FeatureCounts et Limma.
* Puis nous avons intégré cette méthode dans un pipeline afin d'obtenir une entrée compatible avec le fonctionnement de Moomin.
* Nous avons ensuite évalué la qualité des résultats finaux obtenus avec Moomin en les comparant avec ceux obtenus classiquement avec des données RNASeq.
* Ça m'a permis de me familiariser avec les méthodes mises au point par l'équipe et avec des concepts biologiques importants pour la suite.

## Perspectives :

- [x] Valider la méthode avec des données microarrays avec un jeu où l'on a également directement les flux (prévu avec Delphine)
- [x] Intégrer données protéomiques, métabolomiques (Thèse)
- [x] Adapter à un modèle de communauté microbienne (prise en compte de l'environnement)

## Logiciels et bibliothèques:

- Moomin
- Escher (v1.7.0-beta.15)
- COBRA Toolbox (v3.0)
- Cyber-t  (v2.0beta)
- EBSeq (v1.23.1)
- Limma (v3.39.19)
- FeatureCounts (v1.5.3)
- Bowtie2 (v2.3.4.1)
- FastQC (v0.11.5)
- Cutadapt (v1.18)

## Environnement de développement :

- Linux (Ubuntu 18.04)
- R (v3.6)
- Matlab (R2016a)
- Python (v3.7.3)