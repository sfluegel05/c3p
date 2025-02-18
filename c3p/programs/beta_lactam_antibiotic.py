"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: CHEBI:27731 beta-lactam antibiotic
"""
from rdkit import Chem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    A beta-lactam antibiotic contains a beta-lactam ring (a four-membered cyclic amide).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam antibiotic, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for a basic beta-lactam ring
    beta_lactam_pattern = Chem.MolFromSmarts("C1C(=O)N1")

    # Check if basic beta-lactam is present first
    if not mol.HasSubstructMatch(beta_lactam_pattern):
       return False, "Does not contain a beta-lactam ring"

    # Define SMARTS pattern for common beta-lactam ring systems (more general)
    # Penam (Penicillin) ring system (thiazolidine fused to beta-lactam) - Relaxed constraints
    penam_pattern = Chem.MolFromSmarts("S1CC[CH]N2C(C1)C(=O)C2")

    # Cephem (Cephalosporin) ring system (dihydrothiazine fused to beta-lactam) - Relaxed constraints
    cephem_pattern = Chem.MolFromSmarts("S1C=C[CH]N2C(C1)C(=O)C2")
    
    # Carbapenem ring system (5 membered ring fused to beta-lactam) - Relaxed constraints
    carbapenem_pattern = Chem.MolFromSmarts("C1CC[CH]N2C(C1)C(=O)C2")


    # Check for substructure match, starting with most specific ones
    if mol.HasSubstructMatch(penam_pattern):
         return True, "Contains a penam (penicillin) beta-lactam ring system"
    elif mol.HasSubstructMatch(cephem_pattern):
        return True, "Contains a cephem (cephalosporin) beta-lactam ring system"
    elif mol.HasSubstructMatch(carbapenem_pattern):
        return True, "Contains a carbapenem beta-lactam ring system"
    else:
      return True, "Contains a beta-lactam ring in a possible antibiotic system"