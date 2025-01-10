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

    # Check each nitrile group
    for match in nitrile_matches:
        carbon_idx = match[0]  # Index of the carbon in the nitrile group
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        
        # Ensure the carbon is not part of an aromatic ring
        if carbon_atom.GetIsAromatic():
            return False, "Nitrile group is attached to an aromatic carbon"
        
        # Ensure the carbon is part of an aliphatic chain
        if not carbon_atom.GetSymbol() == "C":
            return False, "Nitrile group is not attached to a carbon atom"

        # Check if the carbon is part of an aliphatic chain (not just a single carbon)
        neighbors = carbon_atom.GetNeighbors()
        if len(neighbors) < 2:
            return False, "Nitrile group is not part of an aliphatic chain"
        
        # Check if all neighbors are aliphatic carbons
        for neighbor in neighbors:
            if neighbor.GetSymbol() != "C" or neighbor.GetIsAromatic():
                return False, "Nitrile group is not fully attached to aliphatic carbons"

        # Check for conjugation with aromatic systems
        # Get all atoms within 3 bonds of the nitrile carbon
        env = Chem.FindAtomEnvironmentOfRadiusN(mol, 3, carbon_idx)
        atoms_in_env = set()
        for bond_idx in env:
            bond = mol.GetBondWithIdx(bond_idx)
            atoms_in_env.add(bond.GetBeginAtomIdx())
            atoms_in_env.add(bond.GetEndAtomIdx())
        
        # Check if any atom in the environment is aromatic
        for atom_idx in atoms_in_env:
            if mol.GetAtomWithIdx(atom_idx).GetIsAromatic():
                return False, "Nitrile group is conjugated with aromatic system"

    # Additional checks for molecular complexity
    # Count aromatic rings
    n_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if n_aromatic_rings > 0:
        return False, "Molecule contains aromatic rings"

    # Count rotatable bonds to ensure simple aliphatic structure
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable > 10:
        return False, "Molecule is too complex for simple aliphatic nitrile"

    return True, "Contains a nitrile group attached to an aliphatic carbon chain"