"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: monounsaturated fatty acid
"""
from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.
    A monounsaturated fatty acid (MUFA) is a fatty acid with exactly one double or triple bond
    in the fatty acid chain and singly bonded carbon atoms in the rest of the chain.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group (C(=O)O[H])
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[OH]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"
    
    # Get index of carboxyl carbon
    carboxyl_carbons = mol.GetSubstructMatch(Chem.MolFromSmarts('C(=O)[OH]'))
    if not carboxyl_carbons:
        return False, "No carboxylic acid carbon found"
    carboxyl_carbon_idx = carboxyl_carbons[0]
    
    # Use BFS to traverse the fatty acid chain starting from the carboxyl carbon
    visited = set()
    to_visit = [carboxyl_carbon_idx]
    double_triple_bonds = 0
    
    while to_visit:
        atom_idx = to_visit.pop()
        if atom_idx in visited:
            continue
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        atom_symbol = atom.GetSymbol()
        
        if atom_symbol != 'C':
            # Allow oxygen atoms in carboxyl group and possible hydroxyl groups
            if atom_symbol not in ('O', 'H'):
                return False, f"Non-carbon atom ({atom_symbol}) found in chain"
        
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            bond_type = bond.GetBondType()
            # Count double and triple bonds between carbons
            if bond_type in (Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE):
                if atom.GetAtomicNum() == 6 and nbr.GetAtomicNum() == 6:
                    double_triple_bonds += 1
            to_visit.append(nbr_idx)
    
    if double_triple_bonds != 1:
        return False, f"Expected 1 double or triple bond in chain, found {double_triple_bonds}"
    
    return True, "Molecule is a monounsaturated fatty acid"