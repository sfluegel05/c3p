"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is a polyunsaturated fatty acid with three double bonds and a carboxylic acid group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trienoic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for exactly 3 C=C double bonds in a non-cyclic context
    # Exclude rings explicitly
    double_bond_count = 0
    rings = mol.GetRingInfo().AtomRings()
    
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6):
                # Check if the bond is not part of any ring
                if not any(begin_atom.GetIdx() in ring and end_atom.GetIdx() in ring for ring in rings):
                    double_bond_count += 1
    
    if double_bond_count != 3:
        return False, f"Found {double_bond_count} carbon-carbon double bonds, need exactly 3"

    # Check number of carbon atoms - indicative heuristic for sufficient chain length
    num_carbons = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if num_carbons < 12:
        return False, f"Insufficient number of carbon atoms ({num_carbons}) for a typical trienoic fatty acid"
    
    return True, "Contains carboxylic acid group and three double bonds in a non-cyclic carbon chain"