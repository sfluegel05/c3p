"""
Classifies: CHEBI:23899 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    An icosanoid is a signaling molecule arising from oxidation of C20 essential fatty acids
    such as EPA, AA, and DGLA.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify the longest continuous carbon chain
    def get_longest_carbon_chain_length(mol):
        max_length = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                visited = set()
                stack = [(atom.GetIdx(), 1)]
                while stack:
                    current_atom_idx, length = stack.pop()
                    if length > max_length:
                        max_length = length
                    visited.add(current_atom_idx)
                    current_atom = mol.GetAtomWithIdx(current_atom_idx)
                    for neighbor in current_atom.GetNeighbors():
                        neighbor_idx = neighbor.GetIdx()
                        if neighbor.GetAtomicNum() == 6 and neighbor_idx not in visited:
                            stack.append((neighbor_idx, length + 1))
        return max_length

    longest_chain_length = get_longest_carbon_chain_length(mol)
    if longest_chain_length < 18 or longest_chain_length > 22:
        return False, f"Longest carbon chain length ({longest_chain_length}) not in range for icosanoids (18-22)"
    
    # Check for carboxylic acid group (-COOH)
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for oxidative functional groups
    has_hydroxyl = mol.HasSubstructMatch(Chem.MolFromSmarts("[OX2H]"))  # Hydroxyl group
    has_keto = mol.HasSubstructMatch(Chem.MolFromSmarts("C=O"))         # Keto group
    has_epoxy = mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3]1[OX2][CX3]1"))  # Epoxy group
    
    if not (has_hydroxyl or has_keto or has_epoxy):
        return False, "No oxidative functional groups found (hydroxyl, keto, epoxy)"
    
    # Check for multiple double bonds
    num_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    if num_double_bonds < 2:
        return False, f"Only {num_double_bonds} double bonds found, less than expected for icosanoids"
    
    return True, "Molecule matches key features of icosanoids: C20 backbone, oxidative modifications, multiple double bonds"