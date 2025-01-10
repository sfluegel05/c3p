"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic compound containing two or more amino groups.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Initialize amino group count
    amino_group_count = 0
    
    # SMARTS pattern for amino group (excluding amides, nitro groups, nitriles)
    amino_pattern = Chem.MolFromSmarts("[NX3;H1,H2;!$(N-[!#6])]")  # Nitrogen with 1 or 2 hydrogens, not bonded to non-carbon (e.g., O in amide)
    
    # Find amino groups
    matches = mol.GetSubstructMatches(amino_pattern)
    amino_group_count = len(matches)
    
    # Alternative method: Iterate over atoms to find amino groups
    # for atom in mol.GetAtoms():
    #     if atom.GetAtomicNum() == 7:  # Nitrogen atom
    #         if not atom.IsAromatic():
    #             if atom.GetDegree() <= 3:
    #                 if atom.GetTotalNumHs() >= 1:
    #                     num_bonds = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
    #                     if num_bonds <= 3:
    #                         # Check if not part of amide, nitro, nitrile
    #                         is_amide = any([nbr.GetSymbol() == 'C' and nbr.GetDegree() == 3 and nbr.GetTotalValence() == 4 for nbr in atom.GetNeighbors()])
    #                         if not is_amide:
    #                             amino_group_count += 1
                                
    if amino_group_count >= 2:
        return True, f"Contains {amino_group_count} amino groups"
    else:
        return False, f"Contains {amino_group_count} amino group(s), which is less than 2"