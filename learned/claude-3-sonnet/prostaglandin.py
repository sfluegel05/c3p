"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: CHEBI:35473 prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are naturally occurring compounds derived from the parent C20 acid, prostanoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for cyclopentane ring fused to a cyclopentene ring
    fused_rings_pattern = Chem.MolFromSmarts("[C@H]1C[C@@H]2[C@@]1([H])CC=C2")
    if not mol.HasSubstructMatch(fused_rings_pattern):
        return False, "No fused cyclopentane-cyclopentene ring system found"
    
    # Check for common functional groups
    has_hydroxy = any(atom.GetSymbol() == 'O' and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 1 for atom in mol.GetAtoms())
    has_keto = any(atom.GetSymbol() == 'O' and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 2 for atom in mol.GetAtoms())
    has_carboxyl = any(atom.GetSymbol() == 'O' and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 2 and 
                       any(nbr_atom.GetSymbol() == 'O' and sum(nbr_bond.GetBondTypeAsDouble() for nbr_bond in nbr_atom.GetBonds()) == 1
                           for nbr_atom, nbr_bond in zip(atom.GetNeighbors(), atom.GetBonds()))
                       for atom in mol.GetAtoms())
    
    if not any([has_hydroxy, has_keto, has_carboxyl]):
        return False, "No common prostaglandin functional group found"
    
    # Look for lipid chain (long carbon chain)
    lipid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    lipid_matches = mol.GetSubstructMatches(lipid_pattern)
    if not lipid_matches:
        return False, "No lipid chain found"
    
    # Count rotatable bonds to verify long lipid chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Lipid chain too short"
    
    # Check molecular weight and atom counts
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 400:
        return False, "Molecular weight outside typical range for prostaglandins"
    
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    
    if c_count < 18 or c_count > 25 or o_count < 3 or o_count > 6 or h_count < 25 or h_count > 38:
        return False, "Atom counts outside typical ranges for prostaglandins"
    
    return True, "Contains fused cyclopentane-cyclopentene ring system, common prostaglandin functional groups, and a lipid chain"