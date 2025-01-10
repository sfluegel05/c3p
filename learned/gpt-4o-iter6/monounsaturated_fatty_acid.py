"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid (MUFA) based on its SMILES string.
    MUFAs have one double or triple bond in the fatty acid chain and a carboxylic acid group.

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

    # Define carboxylic acid pattern to confirm presence outside just ring closures
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Define patterns for chain double and triple bonds
    double_bond_pattern = Chem.MolFromSmarts("[#6]=[#6]")  # Represents C=C
    triple_bond_pattern = Chem.MolFromSmarts("[#6]#[#6]")  # Represents C#C

    # Determine non-ring double and triple bonds
    num_double_bonds = len([bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and not (bond.IsInRing())])
    num_triple_bonds = len([bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.TRIPLE and not (bond.IsInRing())])
    
    # Count total unsaturations in chains
    num_unsaturations = num_double_bonds + num_triple_bonds

    # A MUFA should present exactly one instance of unsaturation in chain
    if num_unsaturations != 1:
        return False, f"Found {num_unsaturations} unsaturations in chain, need exactly one"
    
    return True, "Molecule is a monounsaturated fatty acid (one double or triple bond in the chain with carboxylic group)"