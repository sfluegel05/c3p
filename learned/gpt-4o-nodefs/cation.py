"""
Classifies: CHEBI:36916 cation
"""
from rdkit import Chem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a cation, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for atoms with positive formal charges
    is_cationic = False
    reasons = []
    for atom in mol.GetAtoms():
        if atom.GetFormalCharge() > 0:
            is_cationic = True
            reasons.append(f"Atom {atom.GetSymbol()} has positive charge {atom.GetFormalCharge()}")
    
    if is_cationic:
        return True, "Molecule contains positively charged atoms: " + ", ".join(reasons)
    else:
        return False, "No positively charged atoms found in the molecule"

# Example use:
# result = is_cation("C[N+](C)(C)c1ccccc1") # trimethylphenylammonium
# print(result)