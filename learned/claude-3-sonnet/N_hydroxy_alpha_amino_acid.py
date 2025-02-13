"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: N-hydroxy-alpha-amino-acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino acid based on its SMILES string.
    These compounds have a carboxylic acid group and at least one hydrogen attached 
    to the amino group is replaced by a hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-hydroxy-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    # Patterns for different types of N-hydroxy groups
    patterns = [
        # N-hydroxy group (N-OH)
        "[CX4][NX3]([H])[OX2H1]",
        # N,N-dihydroxy group (N(OH)2)
        "[CX4][NX3]([OX2H1])[OX2H1]",
        # Hydroxyimino group (=N-OH)
        "[CX4]~[NX2]~[NX2]~[OX2H1]",
        # Alternative hydroxyimino pattern
        "[CX4]~[NX3]~[CX3]=[NX2][OX2H1]"
    ]
    
    # Check for alpha-amino acid backbone with N-hydroxy group
    found_n_hydroxy = False
    for pattern in patterns:
        substructure = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(substructure):
            found_n_hydroxy = True
            break
            
    if not found_n_hydroxy:
        return False, "No N-hydroxy group found"

    # Verify the structure has an alpha carbon with both carboxyl and N-hydroxy group
    alpha_amino_acid = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1][CX4]~[NX2,NX3]")
    if not mol.HasSubstructMatch(alpha_amino_acid):
        return False, "Not an alpha-amino acid structure"

    # Additional check for correct connectivity
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    for match in carboxyl_matches:
        carboxyl_carbon = mol.GetAtomWithIdx(match[0])
        for neighbor in carboxyl_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # carbon
                # Check if this carbon is connected to nitrogen
                for n_neighbor in neighbor.GetNeighbors():
                    if n_neighbor.GetAtomicNum() == 7:  # nitrogen
                        # Check if nitrogen has oxygen attached
                        for o_neighbor in n_neighbor.GetNeighbors():
                            if o_neighbor.GetAtomicNum() == 8:  # oxygen
                                return True, "Contains N-hydroxy group attached to alpha carbon of amino acid"

    return False, "N-hydroxy group not properly connected to alpha carbon"