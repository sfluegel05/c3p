"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdchem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid (MUFA) based on its SMILES string.
    MUFAs have one double or triple bond in the fatty acid chain and a carboxylic acid group.
    They must have a linear or cyclic configuration.

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

    # Define carboxylic acid pattern to confirm presence
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Define patterns for chain double and triple bonds that are NOT in rings
    double_bond_pattern = Chem.MolFromSmarts("[#6]=[#6]")  # Represents C=C
    triple_bond_pattern = Chem.MolFromSmarts("[#6]#[#6]")  # Represents C#C

    # Detect number of double bonds and ensure they're not counted from carboxylic acid
    num_double_bonds = 0
    num_triple_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == rdchem.BondType.DOUBLE and not bond.IsInRing():
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            # Ensure carbon backbone, ignore carboxylic C=O
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                num_double_bonds += 1
        elif bond.GetBondType() == rdchem.BondType.TRIPLE and not bond.IsInRing():
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            # Ensure carbon backbone
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                num_triple_bonds += 1

    # Count total unsaturations in chain
    num_unsaturations = num_double_bonds + num_triple_bonds

    # A MUFA should present exactly one instance of unsaturation in chain
    if num_unsaturations != 1:
        return False, f"Found {num_unsaturations} unsaturations in chain, need exactly one"
    
    return True, "Molecule is a monounsaturated fatty acid (one non-ring double or triple bond in the chain with a carboxylic group)"