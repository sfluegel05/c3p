"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem

def is_wax(smiles: str):
    """
    Determines if a molecule is classified as a wax based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for ester bond pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester bond present"

    # Check for long carbon chains (at least 16 carbons on either side of ester)
    def has_long_chain(substruct):
        # Count carbon atoms in the substructure
        c_count = sum(1 for atom in substruct.GetAtoms() if atom.GetAtomicNum() == 6)
        return c_count >= 16

    esters = mol.GetSubstructMatches(ester_pattern)
    for ester in esters:
        acyl_chain = Chem.PathToSubmol(mol, ester[:2])
        alkoxy_chain = Chem.PathToSubmol(mol, ester[1:])
        if has_long_chain(acyl_chain) and has_long_chain(alkoxy_chain):
            return True, "Contains an ester bond with long hydrocarbon chains"

    return False, "Ester bond present but no long hydrocarbon chains"