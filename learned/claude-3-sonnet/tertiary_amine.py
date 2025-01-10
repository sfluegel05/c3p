"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: tertiary amine compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule contains a tertiary amine group.
    A tertiary amine has a nitrogen atom with three carbon substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains tertiary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Exclude patterns - functional groups that aren't tertiary amines
    exclude_patterns = [
        Chem.MolFromSmarts('[N]=[*]'),  # Imine, etc
        Chem.MolFromSmarts('[N]#[*]'),  # Nitrile, etc
        Chem.MolFromSmarts('[NX3](=[O,S])[*]'),  # N-oxides
        Chem.MolFromSmarts('[N+](=[O])[O-]'),  # Nitro
        Chem.MolFromSmarts('[NX3][CX3](=[OX1])[#6]'),  # Amide
        Chem.MolFromSmarts('[NX3][CX3](=[SX1])[#6]'),  # Thioamide
        Chem.MolFromSmarts('[NX3][SX4](=[OX1])(=[OX1])[#6]'),  # Sulfonamide
        Chem.MolFromSmarts('[n]'),  # Aromatic nitrogen
        Chem.MolFromSmarts('[N+]'),  # Quaternary nitrogen
    ]

    # Pattern for tertiary amine: N connected to exactly 3 carbons
    tertiary_amine_pattern = Chem.MolFromSmarts('[NX3]([CX4,CH2,CH3,CH])[CX4,CH2,CH3,CH][CX4,CH2,CH3,CH]')
    
    # First check for matches to exclude patterns
    for pattern in exclude_patterns:
        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                # Get the nitrogen atom from the match
                n_idx = match[0]  # Assuming nitrogen is first atom in SMARTS
                # Mark this nitrogen as excluded
                mol.GetAtomWithIdx(n_idx).SetProp('excluded', 'true')

    # Now look for tertiary amine pattern
    if mol.HasSubstructMatch(tertiary_amine_pattern):
        matches = mol.GetSubstructMatches(tertiary_amine_pattern)
        for match in matches:
            n_atom = mol.GetAtomWithIdx(match[0])
            # Skip if this nitrogen was marked as excluded
            if not n_atom.HasProp('excluded'):
                # Check hybridization
                if n_atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                    # Count total bonds to verify it's a tertiary amine
                    if len(n_atom.GetBonds()) == 3:
                        # Count carbon neighbors
                        carbon_neighbors = sum(1 for neighbor in n_atom.GetNeighbors() 
                                            if neighbor.GetAtomicNum() == 6)
                        if carbon_neighbors == 3:
                            return True, "Contains nitrogen atom bonded to exactly three carbon atoms"

    return False, "No tertiary amine group found"