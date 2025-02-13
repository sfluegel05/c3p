"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
from rdkit import Chem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    A beta-lactam antibiotic contains a beta-lactam ring (four-membered cyclic amide), 
    which is a key feature of antibiotics such as penicillins, cephalosporins, 
    monobactams, and carbapenems.

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
    
    # SMARTS pattern for detecting the beta-lactam ring
    beta_lactam_ring_pattern = Chem.MolFromSmarts("C1C(=O)N(C)C1")
    
    if mol.HasSubstructMatch(beta_lactam_ring_pattern):
        # Check for additional structural features typical of beta-lactam antibiotics
        # such as specific side chains or functional groups
        penicillin_pattern = Chem.MolFromSmarts("C1C(=O)N2C(C)C2SC1")
        cephalosporin_pattern = Chem.MolFromSmarts("C1C2=C(C)SC(N2C(N1)=O)C(=O)")
        carbapenem_pattern = Chem.MolFromSmarts("C1C2=C(O)C(=O)N(C2C1)C=O")
        
        if (mol.HasSubstructMatch(penicillin_pattern) or
            mol.HasSubstructMatch(cephalosporin_pattern) or
            mol.HasSubstructMatch(carbapenem_pattern)):
            return True, "Contains beta-lactam ring with known antibiotic structure features"

        return True, "Contains a beta-lactam ring"
    
    return False, "No beta-lactam ring found"