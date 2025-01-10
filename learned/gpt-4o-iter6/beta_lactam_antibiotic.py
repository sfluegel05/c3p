"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
from rdkit import Chem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    A beta-lactam antibiotic typically contains a beta-lactam ring: a 4-membered ring with one nitrogen atom.
    Including penicillins, cephalosporins, and carbapenems.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is identified as a beta-lactam antibiotic, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # Define SMARTS patterns for beta-lactam antibiotics
    beta_lactam_patterns = [
        # General beta-lactam ring (azetidinone)
        Chem.MolFromSmarts("C1CNC(=O)1"), 
        
        # Penicillin core: thiazolidine fused with beta-lactam
        Chem.MolFromSmarts("C1SCCN1C(=O)"),  
        
        # Cephalosporin core: dihydrothiazine fused with beta-lactam
        Chem.MolFromSmarts("C1SCC2=C(N1C(=O))"), 
        
        # Carbapenem core: unsaturated beta-lactam with sulfur
        Chem.MolFromSmarts("C1CC(=O)N1C=O"),  
    ]
    
    # If any pattern matches the molecular structure
    for pattern in beta_lactam_patterns:
        if mol.HasSubstructMatch(pattern):
            return (True, "Contains a beta-lactam ring with antibiotic scaffold")

    return (False, "No beta-lactam ring with known antibiotic scaffold found")