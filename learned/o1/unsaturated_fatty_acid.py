"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
"""
Classifies: unsaturated fatty acid
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid is any fatty acid containing at least one C=C or C#C bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[O;H1]')
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid)
    if not carboxy_matches:
        return False, "No carboxylic acid group found"

    # Consider each carboxylic acid group found
    for match in carboxy_matches:
        carboxyl_carbon_idx = match[0]
        carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)
        
        # Traverse all connected carbon atoms starting from the carboxyl carbon
        visited = set()
        carbon_chain_atoms = []

        def traverse_chain(atom):
            if atom.GetIdx() in visited:
                return
            visited.add(atom.GetIdx())
            if atom.GetAtomicNum() == 6:
                carbon_chain_atoms.append(atom)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 or neighbor.GetAtomicNum() == 1:
                    traverse_chain(neighbor)
                
        # Start traversal from the carboxyl carbon
        traverse_chain(carboxyl_carbon)
        
        # Exclude the carboxyl carbon itself from the chain atoms
        chain_atoms = [atom for atom in carbon_chain_atoms if atom.GetIdx() != carboxyl_carbon_idx]
        
        if not chain_atoms:
            continue  # No carbon chain attached

        # Check for presence of at least one C=C or C#C bond in the chain
        has_double_triple_bond = False
        for bond in mol.GetBonds():
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom in chain_atoms and end_atom in chain_atoms:
                if bond.GetBondType() in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
                    if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                        has_double_triple_bond = True
                        break
        if not has_double_triple_bond:
            continue  # No carbon-carbon double or triple bond in the chain

        # If all checks passed
        return True, "Molecule is an unsaturated fatty acid"

    # If no carboxylic acid group with appropriate chain found
    return False, "Does not meet criteria for unsaturated fatty acid"