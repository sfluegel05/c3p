"""
Classifies: CHEBI:46722 carbonate ester
"""
"""
Classifies: CHEBI:33575 carbonate ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester is a carbonic acid derivative where the hydrogen atoms
    are replaced by organyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carbonate ester pattern (-O-C(=O)-O-) in SMARTS
    carbonate_ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[OX2]")
    matches = mol.GetSubstructMatches(carbonate_ester_pattern)
    
    if not matches:
        return False, "No carbonate ester group found"
    
    # Check that the carbonate ester group is not part of other functional groups
    for match in matches:
        c_atom = mol.GetAtomWithIdx(match[1])
        o1_atom = mol.GetAtomWithIdx(match[0])
        o2_atom = mol.GetAtomWithIdx(match[2])
        
        # Check that carbonyl carbon is not part of amide, ester, or anhydride
        if c_atom.GetTotalNumHs() == 0 and any(nb.GetAtomicNum() == 7 for nb in c_atom.GetNeighbors()):
            continue  # Amide
        if any(nb.GetAtomicNum() == 8 and nb.GetTotalNumHs() == 0 for nb in c_atom.GetNeighbors()):
            continue  # Ester or anhydride
        
        # Check that oxygens are not part of ether, alcohol, or other functional groups
        if any(nb.GetAtomicNum() == 8 for nb in o1_atom.GetNeighbors()) or o1_atom.GetTotalNumHs() > 0:
            continue
        if any(nb.GetAtomicNum() == 8 for nb in o2_atom.GetNeighbors()) or o2_atom.GetTotalNumHs() > 0:
            continue
        
        return True, "Contains carbonate ester group (-O-C(=O)-O-)"
    
    return False, "Identified groups are not carbonate esters"