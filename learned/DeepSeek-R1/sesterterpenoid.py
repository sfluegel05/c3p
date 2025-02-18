"""
Classifies: CHEBI:26660 sesterterpenoid
"""
"""
Classifies: CHEBI_XXXXX sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    Sesterterpenoids are derived from a C25 sesterterpene skeleton, which may be modified.

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
    
    # Check carbon count (allow some modification from base C25)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (23 <= c_count <= 27):
        return False, f"Carbon count {c_count} outside sesterterpenoid range (23-27)"
    
    # Calculate approximate isoprene units (target ~5 for sesterterpenoids)
    num_isoprene = rdMolDescriptors.CalcNumIsopreneUnits(mol)
    if not (4.5 <= num_isoprene <= 5.5):
        return False, f"{num_isoprene:.1f} isoprene units (expected ~5)"
    
    # Check for common terpenoid features (e.g., multiple methyl groups)
    methyl_pattern = Chem.MolFromSmarts("[CH3]")
    methyl_matches = len(mol.GetSubstructMatches(methyl_pattern))
    if methyl_matches < 2:  # Allow some demethylation
        return False, "Insufficient methyl groups for terpenoid structure"
    
    return True, "Meets sesterterpenoid criteria: ~25 carbons, ~5 isoprene units, and methyl groups"