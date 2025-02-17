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
    i.e. a four-membered ring containing one nitrogen and one carbonyl carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a beta-lactam antibiotic, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS for a beta-lactam ring:
    # - [NX2;r4]: a nitrogen in a 4-membered ring with two bonds
    # - [CX4;r4]: a saturated carbon in the ring
    # - [CX3;r4](=O): an sp2 carbon (carbonyl) in the ring with a double bond to oxygen
    beta_lactam_smarts = "[NX2;r4]1[CX4;r4][CX3;r4](=O)[CX4;r4]1"
    beta_lactam_pattern = Chem.MolFromSmarts(beta_lactam_smarts)
    
    if beta_lactam_pattern is None:
        return False, "Failed to generate SMARTS pattern for beta-lactam ring"

    # Check for the presence of the beta-lactam ring substructure
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Molecule contains a beta-lactam ring and is classified as a beta-lactam antibiotic"
    else:
        return False, "No beta-lactam ring found in the molecule"