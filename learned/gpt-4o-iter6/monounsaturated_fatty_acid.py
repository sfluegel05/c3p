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

    # Check for the terminal carboxylic acid group pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No terminal carboxylic acid group found"

    # Calculate number of chain double/triple bonds
    # Single bonds only should exist between carbon atoms in chain, not part of any ring system
    num_unsaturations = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
            start_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            # Ensure neither atom is part of a ring structure
            if not start_atom.IsInRing() and not end_atom.IsInRing():
                num_unsaturations += 1

    if num_unsaturations != 1:
        return False, f"Found {num_unsaturations} unsaturations in chain, need exactly one"
    
    return True, "Molecule is a monounsaturated fatty acid (one double or triple bond in the chain with carboxylic group)"