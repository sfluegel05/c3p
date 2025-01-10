"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: CHEBI:26218 porphyrin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin is a macrocyclic structure with four pyrrole rings connected by methine bridges.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the core porphyrin structure
    # Pattern: 4 pyrrole-like rings connected by methine bridges
    porphyrin_pattern = Chem.MolFromSmarts("c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3")
    if not mol.HasSubstructMatch(porphyrin_pattern):
        return False, "No porphyrin core structure found"

    # Check for at least 4 nitrogen atoms (one in each pyrrole ring)
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if nitrogen_count < 4:
        return False, f"Only {nitrogen_count} nitrogen atoms found, need at least 4"

    # Check ring count - should have at least 5 rings (4 pyrroles + macrocycle)
    ring_info = mol.GetRingInfo()
    if len(ring_info.AtomRings()) < 5:
        return False, "Not enough rings for porphyrin structure"

    # Check molecular weight - porphyrins are typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for porphyrin"

    # Count carbons - porphyrins typically have >20 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons for porphyrin"

    return True, "Contains characteristic porphyrin macrocyclic structure with four pyrrole rings"