// Translations data for Attributor Guide
// This file defines the translations object directly to avoid CORS issues with file:// protocol

const attributorTranslations = {
  "meta": {
    "version": "1.0",
    "languages": ["fr", "en"],
    "default": "fr"
  },
  "nav": {
    "fr": {
      "introduction": "Introduction",
      "prerequisites": "Prérequis",
      "installation": "Installation",
      "interface": "Interface",
      "workflow": "Flux de travail",
      "advanced": "Fonctionnalités avancées",
      "algorithms": "Algorithmes",
      "tutorial": "Tutoriel",
      "troubleshooting": "Dépannage",
      "license": "Licence"
    },
    "en": {
      "introduction": "Introduction",
      "prerequisites": "Prerequisites",
      "installation": "Installation",
      "interface": "Interface",
      "workflow": "Workflow",
      "advanced": "Advanced Features",
      "algorithms": "Algorithms",
      "tutorial": "Tutorial",
      "troubleshooting": "Troubleshooting",
      "license": "License"
    }
  },
  "hero": {
    "fr": {
      "title": "Bienvenue dans Attributor",
      "subtitle": "Un outil puissant développé à l'IPAG (Institut de Planétologie et d'Astrophysique de Grenoble) pour l'analyse et l'attribution des spectres de masse en astronomie et chimie analytique."
    },
    "en": {
      "title": "Welcome to Attributor",
      "subtitle": "A powerful tool developed at IPAG (Institut de Planétologie et d'Astrophysique de Grenoble) for mass spectrometry analysis and molecular attribution in astronomy and analytical chemistry."
    }
  },
  "introduction": {
    "fr": {
      "title": "Introduction",
      "p1": "Attributor est un logiciel spécialisé pour Igor Pro conçu pour aider les chercheurs à :",
      "list": [
        "Construire des formules moléculaires de manière interactive via un tableau périodique",
        "Simuler les distributions isotopiques théoriques avec une précision élevée",
        "Comparer les spectres simulés avec les données expérimentales",
        "Attribuer des formules moléculaires à des pics de masse inconnus",
        "Calibrer les données expérimentales en fonction des formules théoriques",
        "Analyser les écarts en masse et en probabilité"
      ],
      "domain_title": "Domaine d'application",
      "domain_text": "Attributor est particulièrement utile en astrophysique moléculaire, chimie analytique, et spectrométrie de masse haute résolution. Il permet d'identifier et de caractériser des molécules complexes à partir de leurs spectres de masse, en tenant compte des différentes combinaisons isotopiques possibles."
    },
    "en": {
      "title": "Introduction",
      "p1": "Attributor is a specialized software for Igor Pro designed to help researchers:",
      "list": [
        "Interactively build molecular formulas via a periodic table",
        "Simulate theoretical isotopic distributions with high precision",
        "Compare simulated spectra with experimental data",
        "Attribute molecular formulas to unknown mass peaks",
        "Calibrate experimental data based on theoretical formulas",
        "Analyze mass and probability deviations"
      ],
      "domain_title": "Application Domain",
      "domain_text": "Attributor is particularly useful in molecular astrophysics, analytical chemistry, and high-resolution mass spectrometry. It enables the identification and characterization of complex molecules from their mass spectra, taking into account the different possible isotopic combinations."
    }
  },
  "prerequisites": {
    "fr": {
      "title": "Prérequis",
      "steps": {
        "1": {"title": "Igor Pro 64-bit", "desc": "Le logiciel nécessite Igor Pro version 64-bit (WaveMetrics) installé sur votre machine."},
        "2": {"title": "Système d'exploitation", "desc": "Compatible avec Windows 7/10/11 (version 64-bit recommandée)."},
        "3": {"title": "Espace disque", "desc": "Le fichier principal .pxp et ses dépendances nécessitent environ 10-20 Mo d'espace disque."}
      },
      "info_title": "Version d'Igor Pro",
      "info_text": "Attributor a été développé et testé avec les versions récentes d'Igor Pro (8.x et supérieures). Pour une compatibilité optimale, utilisez la dernière version disponible d'Igor Pro 64-bit."
    },
    "en": {
      "title": "Prerequisites",
      "steps": {
        "1": {"title": "Igor Pro 64-bit", "desc": "The software requires Igor Pro 64-bit version (WaveMetrics) installed on your machine."},
        "2": {"title": "Operating System", "desc": "Compatible with Windows 7/10/11 (64-bit version recommended)."},
        "3": {"title": "Disk Space", "desc": "The main .pxp file and its dependencies require approximately 10-20 MB of disk space."}
      },
      "info_title": "Igor Pro Version",
      "info_text": "Attributor was developed and tested with recent versions of Igor Pro (8.x and above). For optimal compatibility, use the latest available version of Igor Pro 64-bit."
    }
  },
  "installation": {
    "fr": {
      "title": "Installation et Lancement",
      "files_title": "Fichiers fournis",
      "procedure_title": "Procédure d'installation",
      "first_launch": "Premier lancement",
      "windows_text": "Lors du premier lancement, les fenêtres suivantes devraient s'ouvrir automatiquement : panel, elaborateur, agregateur, dmvm, AdvancedManager.",
      "warning_title": "Problème de chargement ?",
      "warning_text": "Si les fenêtres ne s'ouvrent pas automatiquement, vérifiez que MainProc.ipf est dans le même dossier que ATTRIBUTOR 2022.pxp et qu'aucune erreur n'est affichée dans la fenêtre History."
    },
    "en": {
      "title": "Installation and Launch",
      "files_title": "Provided Files",
      "procedure_title": "Installation Procedure",
      "first_launch": "First Launch",
      "windows_text": "On first launch, the following windows should open automatically: panel, elaborateur, agregator, dmvm, AdvancedManager.",
      "warning_title": "Loading Problem?",
      "warning_text": "If windows do not open automatically, verify that MainProc.ipf is in the same folder as ATTRIBUTOR 2022.pxp and that no errors are displayed in the History window."
    }
  },
  "interface": {
    "fr": {
      "title": "Interface Utilisateur",
      "screenshot_title": "Description de l'interface (d'après votre capture d'écran)",
      "screenshot_desc": "La capture d'écran montre l'interface complète d'Attributor avec plusieurs fenêtres ouvertes.",
      "main_windows": "Fenêtres principales"
    },
    "en": {
      "title": "User Interface",
      "screenshot_title": "Interface Description (based on your screenshot)",
      "screenshot_desc": "The screenshot shows the complete Attributor interface with multiple windows open.",
      "main_windows": "Main Windows"
    }
  },
  "buttons": {
    "fr": {
      "toggle": "Changer pour l'anglais"
    },
    "en": {
      "toggle": "Switch to French"
    }
  }
};
