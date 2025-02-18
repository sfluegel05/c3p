"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: CHEBI tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    A tetraterpenoid is a terpenoid derived from a tetraterpene.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for tetraterpenoid"

    # Check for isoprene units
    isoprene_pattern = Chem.MolFromSmarts("C=C(C)C=C(C)|C=C(C)C=C(C)C=C(C)")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 4:
        return False, f"Found {len(isoprene_matches)} isoprene units, need at least 4"

    # Check for conjugated double bonds
    conjugated_pattern = Chem.MolFromSmarts("C=C/C=C")
    conjugated_matches = mol.GetSubstructMatches(conjugated_pattern)
    if len(conjugated_matches) < 3:
        return False, f"Found {len(conjugated_matches)} conjugated double bonds, need at least 3"

    # Check for methyl groups
    methyl_pattern = Chem.MolFromSmarts("[CX4H3]")
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    if len(methyl_matches) < 5:
        return False, f"Found {len(methyl_matches)} methyl groups, need at least 5"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 30:
        return False, "Too few carbons for tetraterpenoid"
    if o_count > 15:
        return False, "Too many oxygens for tetraterpenoid"

    return True, "Contains isoprene units, conjugated double bonds, and methyl groups consistent with tetraterpenoid structure"