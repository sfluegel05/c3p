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
    
    # Count carbons in longest chain
    # First find the carboxyl carbon
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"
    
    carboxyl_carbon = carboxyl_matches[0][0]
    
    # Find longest carbon chain starting from carboxyl group
    def count_chain_carbons(atom, visited):
        if atom.GetAtomicNum() != 6:  # Not carbon
            return 0
        max_length = 0
        visited.add(atom.GetIdx())
        for bond in atom.GetBonds():
            next_atom = bond.GetOtherAtom(atom)
            if next_atom.GetIdx() not in visited and next_atom.GetAtomicNum() == 6:
                length = 1 + count_chain_carbons(next_atom, visited)
                max_length = max(max_length, length)
        visited.remove(atom.GetIdx())
        return max_length
    
    chain_length = 1 + count_chain_carbons(mol.GetAtomWithIdx(carboxyl_carbon), set())
    
    if chain_length != 18:
        return False, f"Longest carbon chain has {chain_length} carbons, must be exactly 18"
    
    # Check for rings
    if rdMolDescriptors.CalcNumRings(mol) > 0:
        return False, "Contains rings, must be acyclic"
    
    # Count C=C double bonds in the main chain
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    
    if double_bonds != 2:
        return False, f"Contains {double_bonds} C=C double bonds, must be exactly 2"
    
    # Check that all carbons are either:
    # - part of the main chain
    # - part of -OH substituents
    # - part of -O-CH3 substituents (methoxy groups)
    allowed_substituents = {
        Chem.MolFromSmarts("[OH1]"),  # hydroxyl
        Chem.MolFromSmarts("[O][CH3]"),  # methoxy
        Chem.MolFromSmarts("[O][OH1]"),  # hydroperoxy
    }
    
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            non_h_bonds = len([b for b in atom.GetBonds() if b.GetOtherAtom(atom).GetAtomicNum() != 1])
            if non_h_bonds > 4:  # Maximum possible bonds for carbon
                return False, "Invalid carbon bonding"
    
    # Check for allowed atoms
    allowed_atoms = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Contains atoms other than C, H, and O"
    
    return True, "C18 straight-chain fatty acid with 2 C=C double bonds"