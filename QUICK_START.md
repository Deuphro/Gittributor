# Attributor - Guide Rapide

## 🚀 Installation en 3 étapes

1. **Vérifier les prérequis**
   - Igor Pro 64-bit installé (WaveMetrics)
   - Windows 7/10/11

2. **Préparer les fichiers**
   - Placez `ATTRIBUTOR 2022.pxp` et `MainProc.ipf` dans le même dossier

3. **Lancer**
   - Double-cliquez sur **Igor64.exe**
   - **Fichier → Ouvrir → Ouvrir Experiment...** → Sélectionnez `ATTRIBUTOR 2022.pxp`

---

## 🖥️ Interface Principale

```
┌─────────────────────────────────────────────────────────────┐
│  AdvancedManager (Gestion des données)                           │
├─────────────────┬─────────────────┬─────────────────────────┤
│  dmvm            │  agregateur     │  elaborateur             │
│  (Mass defect)    │  (Comparaison)  │  (Simulation principale) │
├─────────────────┴─────────────────┴─────────────────────────┤
│                                                                 │
│  panel (Tableau périodique + contrôles)                        │
│  ┌─────────────────────────────────────────────────────────┐│
│  │ H  He Li Be ... C  N  O  F  Ne ... (boutons éléments)     ││
│  │                                 Zeromol (Reset)            ││
│  │                                 Add to agregator           ││
│  │                                 Add as factor              ││
│  │                                 Recal / Hard Recal          ││
│  │                                 doMSMS                     ││
│  └─────────────────────────────────────────────────────────┘│
└─────────────────────────────────────────────────────────────┘
```

---

## 🎯 Workflow Standard

### 1️⃣ Charger les données
- Dans **AdvancedManager** → Cliquez sur **Load**
- Sélectionnez votre fichier de spectre de masse

### 2️⃣ Construire la formule moléculaire
- Dans **panel** (tableau périodique)
- Cliquez sur les éléments pour ajouter des atomes
  - Ex: Pour C₆H₁₂O₆ → Cliquez 6× sur **C**, 12× sur **H**, 6× sur **O**
- Ajustez la **charge ionique** si nécessaire (ex: +1 pour [M+H]⁺)

### 3️⃣ Visualiser et comparer
- La simulation apparaît dans **elaborateur**
- Comparez le profil isotopique (rouge/bleu) avec vos données (gris)

### 4️⃣ Calibrer (si nécessaire)
- Cliquez sur **Recal** pour un ajustement automatique
- Ou **Hard Recal** pour un ajustement forcé

### 5️⃣ Ajouter à la comparaison
- Cliquez sur **Add to agregator** pour superposer avec d'autres formules

---

## ⚙️ Contrôles Importants

| Contrôle | Description | Effet |
|----------|-------------|-------|
| **Reset** | Efface la molécule courante | Remet tous les compteurs à 0 |
| **Add to agregator** | Ajoute à la comparaison | Superpose la simulation dans agregateur |
| **Add as factor** | Ajoute comme facteur | Ajoute à la liste des facteurs |
| **Recal** | Calibration automatique | Décale les données expérimentales |
| **Hard Recal** | Calibration forcée | Ajustement plus strict |
| **doMSMS** | Mode MS/MS | Active l'analyse de fragmentation |

### Paramètres de configuration

| Paramètre | Description | Valeur par défaut |
|-----------|-------------|-------------------|
| Isotopic combo probability | Seuil de probabilité pour les combinaisons isotopiques | 0.01 |
| Display relative to maximum | Affichage relatif ou absolu | Désactivé |
| Ion charge | Charge ionique pour le calcul des masses | 0 |

---

## 📊 Fenêtres et leurs Fonctions

| Fenêtre | Purpose | Visualisation |
|---------|---------|--------------|
| **panel** | Interface principale | Tableau périodique, contrôles |
| **elaborateur** | Simulation principale | Spectre simulé vs données |
| **agregateur** | Comparaison multiple | Plusieurs simulations superposées |
| **dmvm** | Analyse du défaut de masse | Mass defect vs m/z |
| **AdvancedManager** | Gestion des données | Chargement des spectres |

---

## 🎓 Exemple Pratique

**Problème**: Vous avez un pic à m/z = 180.0841

**Solution**:
```
1. Dans panel:
   - Cliquez 10× sur C (Carbone)
   - Cliquez 12× sur H (Hydrogène)
   - Cliquez 3× sur O (Oxygène)
   
2. Formule affichée: C10H12O3
   
3. Dans elaborateur:
   - Masse monoisotopique calculée: 180.0841 ✓
   - Profil isotopique: comparez M, M+1, M+2
   
4. Si correspondance:
   - Molécule identifiée !
   
5. Sinon:
   - Essayez d'autres combinaisons
   - Utilisez Recal pour ajuster
```

---

## 🔍 Identification de Molécules

### Conseils rapides

1. **Commencez simple**: Essayez des formules avec peu d'hétéroatomes (O, N, S)
2. **Utilisez dmvm**: Le défaut de masse aide à distinguer les éléments
3. **Comparez les intensités**: Le rapport M/M+1/M+2 est caractéristique
4. **Charge ionique**: N'oubliez pas de régler la charge pour [M+H]⁺, [M-H]⁻, etc.
5. **Seuil de probabilité**: Augmentez-le (ex: 0.05) pour accélérer les calculs

### Exemples de masses courantes

| Formule | Masse monoisotopique | Utilisation typique |
|---------|---------------------|---------------------|
| CH4 | 16.0313 | Méthane |
| H2O | 18.0106 | Eau |
| CO2 | 43.9898 | Dioxyde de carbone |
| C6H12O6 | 180.0634 | Glucose |
| C8H10N4O2 | 194.0804 | Caféine |

---

## ⚠️ Dépannage Rapide

| Problème | Solution |
|----------|----------|
| .pxp ne s'ouvre pas | Vérifiez Igor Pro 64-bit, vérifiez les chemins |
| Fenêtres manquantes | Vérifiez MainProc.ipf dans le même dossier |
| Erreur "wave not found" | Cliquez Reset, puis réessayez |
| Simulation lente | Augmentez le seuil de probabilité à 0.05 |
| Masses incorrectes | Vérifiez la charge ionique, utilisez Recal |
| Interface figée | Attendez 1-2 min, ou réduisez la complexité |

---

## 📚 Documentation Complète

Pour plus de détails, consultez:
- **[Guide complet](attributor-guide.html)** - Documentation détaillée avec tutoriels

---

## 📞 Support

Développeur: François-Régis ORTHOUS-DAUNAY
Affiliation: IPAG (UMR5374), CNRS, UJF Grenoble, France

© 2014-2026 - Licence non exclusive pour usage de recherche uniquement
