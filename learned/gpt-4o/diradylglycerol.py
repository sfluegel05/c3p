"""
Classifies: CHEBI:76578 diradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is a glycerol backbone with two substituent groups 
    (acyl, alkyl, or alk-1-enyl) attached at any two of the three possible positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[O]C[C@@]([O])(C)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Identify and count substituents attached to glycerol backbone
    substituent_pattern = Chem.MolFromSmarts("[CX3](=O)[O][C@H]1[CX2](O)CO1")
    substituent_matches = mol.GetSubstructMatches(substituent_pattern)
    if len(substituent_matches) != 2:
        return False, f"Found {len(substituent_matches)} substituent groups, need exactly 2"
    
    return True, "Contains glycerol backbone with two substituent groups"

# Example test case
smiles_example = "C([C@@](COC(CCCCCC/C=C\C/C=C\C/C=C\CCCCC)=O)(OC(CCCCCCC/C=C\C/C=C\CCCCC)=O)[H])O"
print(is_diradylglycerol(smiles_example))