"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
"""
Classifies: CHEBI:18353 aliphatic nitrile
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile is a molecule containing a nitrile group (C#N) attached to an aliphatic carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the nitrile group (C#N)
    nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)
    if len(nitrile_matches) == 0:
        return False, "No nitrile group (C#N) found"

    # Check if the carbon attached to the nitrile is aliphatic
    for match in nitrile_matches:
        carbon_idx = match[0]  # Index of the carbon in the nitrile group
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        
        # Ensure the carbon is not part of an aromatic ring
        if carbon_atom.GetIsAromatic():
            return False, "Nitrile group is attached to an aromatic carbon"
        
        # Ensure the carbon is part of an aliphatic chain
        if not carbon_atom.GetSymbol() == "C":
            return False, "Nitrile group is not attached to a carbon atom"

    # Check if the molecule is aliphatic (no aromatic rings)
    aromatic_rings = Chem.GetAromaticRings(mol)
    if len(aromatic_rings) > 0:
        return False, "Molecule contains aromatic rings"

    return True, "Contains a nitrile group attached to an aliphatic carbon chain"