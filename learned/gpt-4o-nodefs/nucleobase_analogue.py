"""
Classifies: CHEBI:67142 nucleobase analogue
"""
from rdkit import Chem
from rdkit.Chem import Draw

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expand definition of core purine and pyrimidine structures to include analogues
    nucleobase_analogue_patterns = [
        "c1nncnc1", # Generic pyrimidine
        "c1n[nH]c[nH]c1", # Modified pyrimidine core
        "c1c[nH]c[nH]c1=O", # Pyrimidine with Oxygen
        "n1cnc2c(ncnc12)" # Generic purine
    ]
    
    for pattern in nucleobase_analogue_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return True, f"Structure matches nucleobase analogue pattern: {pattern}"

    # Enhanced functional group recognition may use RDKit molecular properties:
    # - Molecular weight constraints
    # - Number of hydrogen bond donors/acceptors
    # - Specific molecular fingerprints identifying nucleobase-like moieties

    # Finally check for halogen substitutions common in nucleobase analogues
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[F,Cl,Br,I]")):
        return True, "Contains halogen substitution, common in nucleobase analogues"
    
    return False, "Does not match patterns or typical features of nucleobase analogues"