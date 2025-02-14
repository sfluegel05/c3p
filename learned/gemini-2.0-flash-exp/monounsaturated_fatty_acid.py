"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid (MUFA) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a MUFA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O,OH]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Find the carbon chain connected to the carboxyl group
    chain_pattern = Chem.MolFromSmarts("[CX4](=O)[O,OH]-[*;!H0]") # Carbon connected to carboxyl group
    matches = mol.GetSubstructMatches(chain_pattern)

    if not matches:
        return False, "No carbon chain attached to the carboxyl group"

    # Get the first carbon atom connected to the carboxyl group
    start_atom_index = matches[0][0]  # Get index of the carbon in the carboxyl
    
    # Initialize an empty list to store the atoms in the chain
    chain_atoms = [start_atom_index]

    # Recursively follow the chain from the start atom
    def find_chain(atom_index, visited_atoms, chain_atoms):
      atom = mol.GetAtomWithIdx(atom_index)
      neighbors = atom.GetNeighbors()
      for neighbor in neighbors:
        neighbor_index = neighbor.GetIdx()
        if neighbor_index not in visited_atoms:
            visited_atoms.add(neighbor_index)
            bond = mol.GetBondBetweenAtoms(atom_index, neighbor_index)
            if bond.GetBondType() == Chem.BondType.SINGLE and neighbor.GetAtomicNum() == 6:
                  chain_atoms.append(neighbor_index)
                  find_chain(neighbor_index,visited_atoms,chain_atoms)

    visited_atoms = set([start_atom_index])
    find_chain(start_atom_index, visited_atoms, chain_atoms)

    chain_mol = Chem.PathToSubmol(mol, chain_atoms)
    
    # Count double and triple bonds within the chain
    num_double_bonds = len(chain_mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))
    num_triple_bonds = len(chain_mol.GetSubstructMatches(Chem.MolFromSmarts("C#C")))
    total_unsaturations = num_double_bonds + num_triple_bonds

    if total_unsaturations != 1:
        return False, f"Molecule has {total_unsaturations} double/triple bonds in the fatty acid chain, should have exactly 1"

     # Count the number of carbons in the chain
    c_count = len(chain_atoms)

    # Check for fatty acid chain length (at least 4 carbons in the chain, as smallest MUFA is butenoic acid)
    if c_count < 4:
      return False, f"Carbon chain length ({c_count}) is too short for a fatty acid"
    
    return True, "Monounsaturated fatty acid identified"