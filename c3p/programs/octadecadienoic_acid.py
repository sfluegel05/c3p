"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: octadecadienoic acid
Definition: Any straight-chain, C18 polyunsaturated fatty acid having two C=C double bonds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Count total carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count != 18:
        return False, f"Contains {carbon_count} carbons, must be exactly 18"

    # Check for cyclic structures
    if rdMolDescriptors.CalcNumRings(mol) > 0:
        return False, "Contains rings, must be acyclic"
    
    # Look for C18 chain with carboxylic acid and two double bonds
    # This SMARTS pattern matches a C18 chain with carboxylic acid and exactly two C=C double bonds
    # It uses recursive SMARTS to ensure the carbons are only connected to H or chain carbons
    chain_pattern = Chem.MolFromSmarts(
        "[CH3][CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]~"
        "[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]~"
        "[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]~[CH2,CH]C(=O)O"
    )
    
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No C18 straight chain with carboxylic acid found"
    
    # Count C=C double bonds specifically
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    
    if double_bonds != 2:
        return False, f"Contains {double_bonds} C=C double bonds, must be exactly 2"
    
    # Additional check for substituents
    # Count number of non-H bonds for each carbon
    excessive_branching = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            non_h_bonds = sum(1 for bond in atom.GetBonds() 
                            if bond.GetOtherAtom(atom).GetAtomicNum() != 1)
            # Carbons should have max 2 non-H bonds (single bonds)
            # or 3 if part of carboxyl group or double bond
            if non_h_bonds > 3:
                excessive_branching = True
                break
    
    if excessive_branching:
        return False, "Contains branching - must be straight chain"
    
    # Check that we don't have other functional groups besides -OH
    allowed_atoms = set([1, 6, 8])  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, "Contains atoms other than C, H, and O"
    
    return True, "C18 straight-chain fatty acid with 2 C=C double bonds"