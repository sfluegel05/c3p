"""
Classifies: CHEBI:67142 nucleobase analogue
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define more comprehensive nucleobase analogue patterns
    nucleobase_analogue_patterns = [
        "c1[nH]c[nH]c[nH]n1",  # Common pyrimidine modifications
        "n1ccnc2[nH]c[nH]c12",  # Common purine modifications
        "c1nc(=O)[nH]cc1",  # Cytosine derivatives
        "c1nc[nH]c(=O)n1",  # Uracil derivatives
        "n1cc[nH]c2[nH]cnc12",  # Imidazole ring extensions
        "n1cnc2[nH]cnc12"  # Modified purine base
    ]

    for pattern in nucleobase_analogue_patterns:
        substruct = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(substruct):
            return True, f"Structure matches nucleobase analogue pattern: {pattern}"
    
    # Check for specific functional groups or chemical features
    rotatable_bonds = AllChem.CalcNumRotatableBonds(mol)
    heavy_atoms = Chem.rdMolDescriptors.CalcNumHeavyAtoms(mol)
    if 2 <= rotatable_bonds <= 5 and 10 <= heavy_atoms <= 30:
        return False, "Does not match specific structural patterns of nucleobase analogues, but typical physical properties align"

    return False, "Does not match patterns or typical features of nucleobase analogues"