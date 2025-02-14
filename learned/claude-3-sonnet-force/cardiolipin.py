"""
Classifies: CHEBI:28494 cardiolipin
"""
"""
Classifies: CHEBI:36103 cardiolipin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    A cardiolipin is a phosphatidylglycerol composed of two molecules of phosphatidic acid
    covalently linked to a molecule of glycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiolipin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone with 3 oxygen atoms
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for 2 phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("[P+](=O)([O-])[O-]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) != 2:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need exactly 2"

    # Check for 4 fatty acid chains of appropriate length and unsaturation
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 4:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify chain lengths
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    # Check for double bonds in fatty acid chains
    has_double_bonds = any(bond.GetIsAromatic() and bond.GetBondType() == Chem.BondType.DOUBLE for bond in mol.GetBonds())
    if not has_double_bonds:
        return False, "No double bonds found in fatty acid chains"

    # Check molecular weight and composition
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000 or mol_wt > 2000:
        return False, "Molecular weight outside the expected range for cardiolipins"

    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)

    if c_count < 50 or o_count < 15 or p_count != 4:
        return False, "Composition inconsistent with cardiolipins"

    return True, "Contains glycerol backbone with 2 phosphatidic acid moieties attached"