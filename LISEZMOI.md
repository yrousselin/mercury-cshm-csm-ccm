# Analyse de la géométrie de coordination (CShM, CSM, CCM) pour Mercury

*[English version: README.md](README.md)*

Script Python pour **Mercury** (CCDC). Il décrit la forme, la symétrie et la chiralité de la **sphère de coordination des métaux** de la structure ouverte. Les calculs sont faits par la bibliothèque **cosymlib**, et les résultats s'affichent sous forme de rapport dans Mercury.

## Ce que fait le script

Pour chaque centre métallique ayant de 2 à 12 atomes donneurs, le rapport donne :

- **Continuous Shape Measures (CShM) :** l'écart entre le polyèdre de coordination réel et chaque polyèdre idéal de référence pour ce nombre de coordination (octaèdre, prisme trigonal, antiprisme carré…). La référence de plus faible valeur, surlignée en vert, est la géométrie la plus proche. **CShM = 0** correspond à une forme parfaite.
- **Continuous Symmetry Measures (CSM) :** l'écart de la sphère de coordination à des groupes ponctuels comme C₂, C₃, Cₛ ou Cᵢ. **CSM = 0** correspond à une symétrie parfaite.
- **Continuous Chirality Measure (CCM) :** la plus petite CSM parmi les groupes impropres (Cₛ, Cᵢ, S₄, S₈, S₁₂). **CCM = 0** signifie achiral ; **CCM > 0** signifie chiral. Elle est donnée deux fois :
  - *polyèdre seul* : tous les atomes donneurs sont traités comme identiques (chiralité de la forme) ;
  - *avec la nature des atomes* : N, O, etc. sont traités comme des atomes différents.

**Quels atomes sont pris en compte**

- Si vous **sélectionnez des atomes** dans Mercury avant de lancer le script (exactement un métal et ses atomes donneurs), seul ce site est analysé.
- Sinon, chaque métal indépendant par symétrie est analysé. Sa sphère de coordination complète est **reconstruite à partir de la maille et des opérateurs de symétrie** : le script fonctionne donc aussi pour les **MOF, les polymères de coordination et les métaux sur position spéciale** (par exemple Zn dans le ZIF-8, où l'unité asymétrique ne contient qu'un des quatre N). Les atomes donneurs sont les atomes non hydrogène et non métalliques plus proches que la somme des rayons covalents + 0,45 Å. Si la structure n'a pas de maille, ce sont les liaisons définies dans Mercury qui sont utilisées.
- Le rapport donne, pour chaque métal, la liste des atomes donneurs et leurs distances, pour vérifier la sphère de coordination d'un coup d'œil.

## Prérequis

| Quoi | Pourquoi | Où le trouver |
|---|---|---|
| **Mercury avec la CSD Python API** | exécute le script | [CCDC](https://www.ccdc.cam.ac.uk/) (licence CSD nécessaire) |
| **cosymlib** | calcule les CShM, CSM et CCM | [documentation de cosymlib](https://cosymlib.readthedocs.io/en/latest/) |
| numpy | géométrie | inclus dans la CSD Python API |

## Installation

### Étape 1 : installer cosymlib

Le script a besoin d'un Python dans lequel `import cosymlib` fonctionne. Quand cosymlib n'est pas disponible dans Mercury, le script appelle le Python indiqué par `PYTHON_EXTERNE` (étape 2). Deux possibilités :

#### Option A (la plus simple) : installer cosymlib dans le Python de la CSD

Appelez pip **par le chemin complet** du Python de la CSD : cosymlib s'installe alors dans ce Python, quel que soit l'environnement conda actif dans votre terminal.

```
/Users/<nom>/CCDC/ccdc-software/csd-python-api/miniconda/bin/python -m pip install cosymlib
/Users/<nom>/CCDC/ccdc-software/csd-python-api/miniconda/bin/python -c "import cosymlib; print('cosymlib OK')"
```

(Windows : `"C:\Program Files\CCDC\ccdc-software\csd-python-api\miniconda\python.exe" -m pip install cosymlib`, en adaptant le chemin à votre installation.)

cosymlib existe en version précompilée pour Python 3.10 à 3.12 sous Linux, Windows et Mac Apple Silicon. Sur un Mac Intel, pip doit compiler cosymlib et a besoin d'un compilateur Fortran.

> ⚠️ **Une mise à jour de la CSD réinstalle ce Python et efface cosymlib.** Après chaque mise à jour, relancez la commande `pip install`.

#### Option B : installer cosymlib dans un environnement à part

Comme le recommande la [page d'installation de cosymlib](https://cosymlib.readthedocs.io/en/latest/install.html) ; cet environnement résiste aux mises à jour de la CSD.

1. Installez **Miniconda** ou **Anaconda** si ce n'est pas déjà fait.
2. Ouvrez un terminal :
   - macOS : l'application **Terminal** ;
   - Windows : **Anaconda Prompt**, depuis le menu Démarrer.
3. Tapez ces commandes une par une, en appuyant sur Entrée après chacune :

   ```
   conda create -n cosym_env python=3.11 -y
   conda activate cosym_env
   pip install cosymlib
   python -c "import cosymlib; print('cosymlib OK')"
   ```

   La dernière ligne doit afficher `cosymlib OK`. En cas d'échec, consultez les remarques par système sur la page d'installation de cosymlib (compilateurs sous Windows, Mac Apple Silicon).
4. Toujours dans ce terminal, notez le chemin complet de ce Python :
   - macOS : `which python`, qui donne par exemple `/Users/<nom>/miniconda3/envs/cosym_env/bin/python`
   - Windows : `where python`, qui donne par exemple `C:\Users\<nom>\miniconda3\envs\cosym_env\python.exe`

### Étape 2 : indiquer ce Python au script

Ouvrez `CShM_CSM_CCM_analysis.py` dans un éditeur de texte (TextEdit en mode texte brut, Bloc-notes, ou un éditeur de code). Remplacez le chemin de la ligne `PYTHON_EXTERNE = ...` (**ligne 10**) par le Python de l'étape 1 :

```python
PYTHON_EXTERNE = "/Users/<nom>/CCDC/ccdc-software/csd-python-api/miniconda/bin/python"   # option A, macOS
PYTHON_EXTERNE = "/Users/<nom>/miniconda3/envs/cosym_env/bin/python"        # option B, macOS
PYTHON_EXTERNE = "C:/Users/<nom>/miniconda3/envs/cosym_env/python.exe"      # option B, Windows
```

Sous Windows, écrivez le chemin avec des barres obliques `/`, comme ci-dessus.

### Étape 3 : ajouter le script à Mercury

1. Placez le script dans un dossier réservé à vos scripts Mercury.
2. Dans Mercury, ouvrez le menu **CSD Python API** → **Options…** → **Add Location**, puis choisissez ce dossier.
3. Le script apparaît alors dans le menu **CSD Python API**.

## Utilisation

1. Ouvrez une structure contenant au moins un métal.
2. Facultatif : sélectionnez un métal et ses atomes donneurs.
3. Lancez **CSD Python API → CShM_CSM_CCM_analysis**.

## Limites

- Nombres de coordination de 2 à 12 uniquement. Les autres sites sont signalés comme ignorés.
- Les atomes donneurs sont trouvés par un critère de distance (somme des rayons covalents + `BOND_TOL` = 0,45 Å, au plus `R_SEARCH` = 3,3 Å). Pour des contacts métal–ligand longs ou ambigus (carboxylates chélatants, contacts axiaux faibles), ajustez ces deux valeurs en tête du script, ou sélectionnez vous-même les atomes donneurs. Les contacts métal–métal ne sont pas comptés.
- Pour les structures désordonnées, vérifiez la liste des donneurs dans le rapport : des atomes de plusieurs composantes de désordre peuvent être comptés.

## Dépannage

Le rapport affiche le message d'erreur complet. Les cas les plus fréquents :

| Message dans le rapport | Cause | Solution |
|---|---|---|
| `Zn1 : CN 1 ignoré` (ou nombre de coordination trop faible) | ancienne version du script, qui ne lisait que l'unité asymétrique | utiliser la version actuelle |
| `ModuleNotFoundError: No module named 'cosymlib'` | cosymlib n'est pas installé dans le Python de `PYTHON_EXTERNE` (souvent après une mise à jour de la CSD) | lancer la commande `pip install` indiquée à la fin du message |
| `processus tué par le signal 11` | cosymlib a planté sur une mesure | les autres mesures sont quand même données ; celle qui a échoué est signalée en rouge |

## Crédits et références

- **cosymlib**, développé par le groupe Electronic Structure and Symmetry (Universitat de Barcelona), fait tous les calculs de forme, de symétrie et de chiralité : <https://cosymlib.readthedocs.io/>. Merci de citer cosymlib et les articles d'origine des mesures continues de forme et de symétrie si vous publiez des résultats.
- Le script utilise la CSD Python API du Cambridge Crystallographic Data Centre (CCDC).
- Écrit par Yoann Rousselin avec l'aide de Claude (Anthropic) comme assistant de programmation. Ce script fait le lien entre Mercury et cosymlib ; il ne calcule pas lui-même les mesures.

## Licence

MIT. Voir [LICENSE](LICENSE). cosymlib est distribuée sous sa propre licence.

## Auteur

**Yoann Rousselin**, ingénieur de recherche CNRS (cristallographie), ICMUB, Dijon.
Questions ou signalement d'un problème : ouvrez une *Issue* sur cette page.
