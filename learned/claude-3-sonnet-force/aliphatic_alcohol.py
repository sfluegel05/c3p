"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: CHEBI:16991 aliphatic alcohol
An alcohol derived from an aliphatic compound.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OX1H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Check for aliphatic carbon chain
    aliphatic_pattern = Chem.MolFromSmarts("[CH4]")
    if not mol.HasSubstructMatch(aliphatic_pattern):
        return False, "No aliphatic carbon chain found"

    # Check for aromatic rings
    if mol.GetNumAtomRings() > 0:
        aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() if Chem.rdchem.IsAromaticRing(mol, ring)]
        if len(aromatic_rings) > 0:
            return False, "Molecule contains aromatic rings, which is not allowed for aliphatic alcohols"

    return True, "Molecule contains a hydroxyl group and an aliphatic carbon chain"