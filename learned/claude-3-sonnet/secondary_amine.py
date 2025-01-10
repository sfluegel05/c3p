"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: CHEBI:32854 secondary amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule contains a secondary amine group.
    A secondary amine has exactly one N-H bond and two carbon atoms attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a secondary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to the molecule
    mol = Chem.AddHs(mol)

    # Find all nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen atom
            # Get number of hydrogens
            n_hydrogens = atom.GetTotalNumHs()
            
            # Get number of carbons directly attached
            carbon_neighbors = sum(1 for neighbor in atom.GetNeighbors() 
                                if neighbor.GetAtomicNum() == 6)
            
            # Get total number of non-hydrogen neighbors
            heavy_neighbors = len([n for n in atom.GetNeighbors()])
            
            # Check if nitrogen is part of special groups we want to exclude
            in_amide = mol.HasSubstructMatch(Chem.MolFromSmarts("[NX3]C(=O)"))
            in_nitro = mol.HasSubstructMatch(Chem.MolFromSmarts("[N+](=O)[O-]"))
            in_nitroso = mol.HasSubstructMatch(Chem.MolFromSmarts("[N]=[O]"))
            in_imine = mol.HasSubstructMatch(Chem.MolFromSmarts("[NX2]=[C,N]"))
            in_azide = mol.HasSubstructMatch(Chem.MolFromSmarts("[N-][N+]#N"))
            
            # Conditions for secondary amine:
            # 1. Exactly one hydrogen
            # 2. Exactly two carbons attached
            # 3. Total of 3 neighbors (2C + 1H)
            # 4. Not part of special groups
            if (n_hydrogens == 1 and 
                carbon_neighbors == 2 and 
                heavy_neighbors == 2 and 
                not any([in_amide, in_nitro, in_nitroso, in_imine, in_azide])):
                return True, "Contains secondary amine group (N-H with two carbons attached)"
            
            # Special case: if we have 1 hydrogen and 2 total heavy neighbors,
            # and both are carbons, it's also a secondary amine
            if (n_hydrogens == 1 and 
                heavy_neighbors == 2 and 
                carbon_neighbors == 2 and 
                not any([in_amide, in_nitro, in_nitroso, in_imine, in_azide])):
                return True, "Contains secondary amine group (N-H with two carbons attached)"

    return False, "No secondary amine group found"