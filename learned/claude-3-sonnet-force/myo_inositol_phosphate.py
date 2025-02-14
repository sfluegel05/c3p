"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: CHEBI:27767 myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol phosphate in which the inositol component has myo-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for inositol backbone pattern (any inositol ring)
    inositol_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@H]([C@H]([C@H]([C@H]1O)O)O)O)O)O")
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_matches:
        return False, "No inositol backbone found"

    # Look for phosphate groups attached to the inositol ring
    inositol_phosphate_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@H]([C@H]([C@H]([C@H]1OP(=[O])([O-])[O-])O)O)O)O")
    inositol_phosphate_matches = mol.GetSubstructMatches(inositol_phosphate_pattern)
    if not inositol_phosphate_matches:
        return False, "No phosphate groups attached to the inositol ring"

    # Count number of phosphate groups
    n_phosphates = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15 and atom.GetFormalCharge() == 0)

    # Check for at least two phosphate groups
    if n_phosphates < 2:
        return False, "Requires at least two phosphate groups"

    # Exclude glycerophospholipids and glycosphingolipids
    glycerol_pattern = Chem.MolFromSmarts("C(C)(C)(OC)")
    sphingosine_pattern = Chem.MolFromSmarts("C(N)(C)(C)")
    if mol.HasSubstructMatch(glycerol_pattern) or mol.HasSubstructMatch(sphingosine_pattern):
        return False, "Molecule is a glycerophospholipid or glycosphingolipid"

    # Exclude molecules with fatty acid chains
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "Molecule contains fatty acid chains"

    return True, f"Contains myo-inositol backbone with {n_phosphates} phosphate group(s)"