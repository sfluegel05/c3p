"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: CHEBI:17741 volatile organic compound
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on its SMILES string.
    A VOC is defined as any organic compound with an initial boiling point ≤ 250°C at 101.3 kPa.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a VOC, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is organic (contains carbon)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count == 0:
        return False, "Not an organic compound (no carbon atoms)"

    # Calculate molecular weight
    mol_wt = Descriptors.MolWt(mol)

    # Check for molecular weight as a proxy for boiling point
    # Adjusted threshold to be more lenient
    if mol_wt > 400:  # Adjusted threshold based on typical VOC molecular weights
        return False, f"Molecular weight ({mol_wt:.2f}) too high for VOC"

    # Check for heavy atoms (non-hydrogen atoms)
    heavy_atom_count = mol.GetNumHeavyAtoms()
    if heavy_atom_count > 30:  # Adjusted threshold for VOCs
        return False, f"Too many heavy atoms ({heavy_atom_count}) for VOC"

    # Check for functional groups that might increase boiling point
    # (e.g., multiple hydroxyl groups, long carbon chains)
    hydroxyl_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2H]")))
    if hydroxyl_count > 3:  # Relaxed restriction
        return False, f"Too many hydroxyl groups ({hydroxyl_count}) for VOC"

    # Check for long carbon chains
    long_chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    if mol.HasSubstructMatch(long_chain_pattern):
        return False, "Long carbon chain detected, likely not a VOC"

    # If all checks pass, classify as VOC
    return True, f"Molecular weight ({mol_wt:.2f}) and structure consistent with VOC criteria"