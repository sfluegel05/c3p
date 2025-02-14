"""
Classifies: CHEBI:167164 mineral nutrient
"""
"""
Classifies: CHEBI:33677 Mineral nutrient
A mineral that is an inorganic nutrient which must be ingested and absorbed in adequate amounts to satisfy 
a wide range of essential metabolic and/or structural functions in the human body.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

# List of common mineral nutrient elements (not exhaustive)
mineral_nutrient_elements = ['K', 'Ca', 'Mg', 'Na', 'Ba', 'Zn', 'Fe', 'Cs', 'Al', 'Sb', 'Pd', 'La']

def is_mineral_nutrient(smiles: str):
    """
    Determines if a molecule is a mineral nutrient based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mineral nutrient, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if the molecule is inorganic
    if Descriptors.MolWt(mol) < 200:
        # Inorganic molecules are typically small (<200 Da)
        
        # Check for presence of mineral nutrient elements
        elements = set([atom.GetSymbol() for atom in mol.GetAtoms()])
        if any(elem in mineral_nutrient_elements for elem in elements):
            
            # Check for common mineral nutrient counterions
            counterions = [item for item in smiles.split('.') if len(item) <= 3 and any(char.isalpha() for char in item)]
            common_counterions = ['[Cl-]', '[O-]', '[OH-]', '[H+]', '[O-]S([O-])(=O)=O', '[O-]C([O-])=O', '[O-][N+]([O-])=O']
            if any(counterion in common_counterions for counterion in counterions):
                return True, "Contains mineral nutrient element(s) and common counterion(s)"
            else:
                return False, "Contains mineral nutrient element(s) but uncommon counterion(s)"
        else:
            return False, "Does not contain mineral nutrient element(s)"
    else:
        return False, "Molecular weight too high for an inorganic mineral nutrient"