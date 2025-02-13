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
        Chem.MolFromSmarts('[NX3](=[OX1])[#6]'),  # N-oxide
        Chem.MolFromSmarts('[NX3](=[SX1])[#6]'),  # N-thiooxide
        Chem.MolFromSmarts('[NX3][CX3](=[OX1])[#6]'),  # Amide
        Chem.MolFromSmarts('[NX3][CX3](=[SX1])[#6]'),  # Thioamide
        Chem.MolFromSmarts('[NX3][SX4](=[OX1])(=[OX1])[#6]'),  # Sulfonamide
        Chem.MolFromSmarts('[n]'),  # Aromatic nitrogen
        Chem.MolFromSmarts('[N+]'),  # Quaternary nitrogen
        Chem.MolFromSmarts('[N-]'),  # Negatively charged nitrogen
        Chem.MolFromSmarts('[NX3]=[NX2]'),  # Diazo
        Chem.MolFromSmarts('[NX3]=[CX3]'),  # Imine
        Chem.MolFromSmarts('[NX3]#[CX2]'),  # Nitrile
        Chem.MolFromSmarts('[N]1[C](=O)[C]1'),  # Aziridine
        Chem.MolFromSmarts('[N]1[C](=O)[CH2][C]1'),  # β-lactam
        Chem.MolFromSmarts('[N+](=[O])[O-]'),  # Nitro
    ]

    # Mark atoms to exclude
    for pattern in exclude_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                n_idx = match[0]  # Assuming nitrogen is first atom in SMARTS
                mol.GetAtomWithIdx(n_idx).SetProp('excluded', 'true')

    # Find all nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen
            # Skip if marked as excluded
            if atom.HasProp('excluded'):
                continue
                
            # Skip if formal charge != 0
            if atom.GetFormalCharge() != 0:
                continue

            # Count carbon neighbors
            carbon_neighbors = 0
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon
                    carbon_neighbors += 1

            # Check if it's a tertiary amine
            if carbon_neighbors == 3:
                # Verify total number of bonds is 3
                if len(atom.GetBonds()) == 3:
                    return True, "Contains a nitrogen atom bonded to exactly three carbons"

    return False, "No tertiary amine group found"