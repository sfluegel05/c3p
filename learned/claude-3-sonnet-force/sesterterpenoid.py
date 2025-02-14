"""
Classifies: CHEBI:26660 sesterterpenoid
"""
"""
Classifies: CHEBI:52348 Sesterterpenoids
A sesterterpenoid is any terpenoid derived from a sesterterpene (C25 backbone). 
The term includes compounds in which the C25 skeleton has been rearranged or modified by the removal of atoms (generally methyl groups).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for terpenoid-like structure
    n_rings = mol.GetRingInfo().NumRings()
    if n_rings < 3:
        return False, "Lacks terpenoid-like ring system"
    
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_carbons < 20 or n_carbons > 30:
        return False, "Number of carbons outside expected range for sesterterpenoid"
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 600:
        return False, "Molecular weight outside expected range for sesterterpenoid"
    
    # Look for sesterterpene-like backbone patterns
    sesterterpene_patterns = [
        "[C&Ring1;R1][C&Ring1;R1][C&Ring1;R1]CC[C&Ring1;R1][C&Ring1;R1][C&Ring1;R1]CC[C&Ring1;R1][C&Ring1;R1][C&Ring1;R1]CC[C&Ring1;R1][C&Ring1;R1][C&Ring1;R1]CC[C&Ring1;R1][C&Ring1;R1][C&Ring1;R1]",
        "[C&Ring1;R1][C&Ring1;R1][C&Ring1;R1]CC[C&Ring1;R1][C&Ring1;R1][C&Ring1;R1]CC[C&Ring1;R1][C&Ring1;R1][C&Ring1;R1]CC[C&Ring1;R1][C&Ring1;R1][C&Ring1;R1]CC[C&Ring1;R1][C&Ring1;R1][C&Ring1;R1]C",
        "[C&Ring1;R1][C&Ring1;R1][C&Ring1;R1]CC[C&Ring1;R1][C&Ring1;R1][C&Ring1;R1]CC[C&Ring1;R1][C&Ring1;R1][C&Ring1;R1]CC[C&Ring1;R1][C&Ring1;R1][C&Ring1;R1]CC[C&Ring1;R1][C&Ring1;R1][C&Ring1;R1]O"
    ]
    
    backbone_match = False
    for pattern in sesterterpene_patterns:
        sesterterpene_backbone = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(sesterterpene_backbone):
            backbone_match = True
            break
    
    if not backbone_match:
        return False, "No sesterterpene-like backbone found"
    
    return True, "Contains terpenoid-like ring system, appropriate number of carbons and molecular weight, and sesterterpene-like backbone"