"""
Classifies: CHEBI:35359 carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    A carboxamidine is defined as having the structure RC(=NR)NR2,
    accounting for variability in R groups and structural contexts.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved SMARTS pattern for carboxamidine
    # This pattern checks for a central carbon double-bonded to nitrogen,
    # single-bonded to another nitrogen which can further bond to any R group including R, H, or alkyl.
    carboxamidine_pattern = Chem.MolFromSmarts('[CX3](=[NX2])[NX3][NX3]')
    
    # Check for matches
    if mol.HasSubstructMatch(carboxamidine_pattern):
        return True, "Contains carboxamidine (RC(=NR)NR2) group"
    
    return False, "No carboxamidine structure found"

# Testing with examples
print(is_carboxamidine("CC(Oc1c(Cl)cccc1Cl)C1=NCCN1"))  # Expected: True for lofexidine