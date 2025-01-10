"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
"""
Classifies: L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Multiple SMARTS patterns to catch different representations of alpha-amino acids
    # Including charged and neutral forms
    patterns = [
        # Standard neutral form
        '[NX3;H2][CX4;H1][CX3](=[OX1])[OX2H1]',
        # Charged form (NH3+ and COO-)
        '[NX4;H3+][CX4;H1][CX3](=[OX1])[OX1-]',
        # Mixed forms
        '[NX3;H2][CX4;H1][CX3](=[OX1])[OX1-]',
        '[NX4;H3+][CX4;H1][CX3](=[OX1])[OX2H1]'
    ]

    found_match = False
    for pattern in patterns:
        aa_pattern = Chem.MolFromSmarts(pattern)
        matches = mol.GetSubstructMatches(aa_pattern)
        
        for match in matches:
            # Get the alpha carbon (second atom in pattern)
            alpha_carbon = mol.GetAtomWithIdx(match[1])
            
            # Check if alpha carbon is chiral
            if alpha_carbon.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
                continue
                
            # Check for proper tetrahedral geometry
            neighbors = [n.GetAtomicNum() for n in alpha_carbon.GetNeighbors()]
            if len(set(neighbors)) < 3:  # Need at least 3 different atoms for chirality
                continue
            
            # Get canonical SMILES with stereochemistry
            atom_smiles = Chem.MolToSmiles(mol, rootedAtAtom=match[1])
            
            # For L-amino acids, we expect @@ in most cases
            # This works because in the canonical SMILES, RDKit will consistently
            # represent L-amino acids with @@ chirality
            if '@@' in atom_smiles:
                found_match = True
                return True, "Found L-alpha-amino acid center with correct stereochemistry"
    
    if not found_match:
        # Check if we found any amino acid pattern at all
        for pattern in patterns:
            aa_pattern = Chem.MolFromSmarts(pattern)
            if mol.HasSubstructMatch(aa_pattern):
                return False, "Found alpha-amino acid pattern but incorrect or unspecified stereochemistry"
        
        return False, "No alpha-amino acid pattern found"

    return False, "No L-alpha-amino acid stereocenter found"