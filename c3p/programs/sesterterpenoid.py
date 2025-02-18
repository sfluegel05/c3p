"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    A sesterterpenoid is derived from a sesterterpene, often featuring rearranged or modified C25 backbones,
    with additional allowance for functional diversity and structural complexity.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Allow a wider range of carbon atoms to capture the diversity
    # Sesterterpenoids are derived from C25 hydrocarbon chains, so we allow up to 67 to accommodate extensive modifications
    if c_count < 20 or c_count > 67:
        return False, f"Carbon count of {c_count} is outside typical range for sesterterpenoids"

    # Check for terpenoid-like patterns and multiple rings
    # Look for patterns of isoprene units (C5H8)
    isoprene_pattern = Chem.MolFromSmarts("C[CH]([CH3])[CH]=[CH2]")
    ring_count = mol.GetRingInfo().NumRings()
    
    if not mol.HasSubstructMatch(isoprene_pattern) or ring_count < 1:
        return False, "Lacks terpenoid-like isoprene units or insufficient structural rings to qualify as a sesterterpenoid"

    # Assess functional group diversity loosely; allow for alcohols, ketones, esters, etc.
    # Rather than count specific functions, check the presence of heteroatoms indicating functional diversity
    heteroatoms_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [7, 8, 16])
    if heteroatoms_count < 1:
        return False, "Lack of heteroatoms indicating insufficient functional diversity"

    return True, "Contains characteristics typical of a sesterterpenoid with a modified and flexible outlook on carbon skeleton and functional diversity"