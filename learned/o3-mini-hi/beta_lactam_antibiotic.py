"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: Beta-lactam antibiotics
Definition: An organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.
The function is_beta_lactam_antibiotic checks if the given SMILES string contains a beta-lactam ring.
"""

from rdkit import Chem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    We classify a molecule as a beta-lactam antibiotic if it contains a beta-lactam ring,
    i.e. a 4-membered ring containing one nitrogen and one carbonyl carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains at least one beta-lactam ring, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define an improved SMARTS pattern for a beta-lactam ring:
    # The pattern requires:
    #   - a four-membered ring (R4),
    #   - one nitrogen (atomic number 7),
    #   - three carbons (atomic number 6), where one of the carbons is double-bonded to oxygen.
    # This SMARTS is flexible to allow stereochemical annotations and common representations.
    beta_lactam_smarts = "[#7;R4]1[#6;R4][#6;R4](=O)[#6;R4]1"
    beta_lactam_pattern = Chem.MolFromSmarts(beta_lactam_smarts)
    
    if beta_lactam_pattern is None:
        return False, "Failed to generate SMARTS pattern for beta-lactam ring"
    
    # Check for the presence of the beta-lactam ring substructure.
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Molecule contains a beta-lactam ring and is classified as a beta-lactam antibiotic"
    else:
        return False, "No beta-lactam ring found in the molecule"