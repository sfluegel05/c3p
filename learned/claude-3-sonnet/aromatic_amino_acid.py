"""
Classifies: CHEBI:33856 aromatic amino acid
"""
"""
Classifies: aromatic amino acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid must have an aromatic ring and an amino acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aromatic ring
    has_aromatic = False
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            has_aromatic = True
            break
            
    if not has_aromatic:
        return False, "No aromatic ring found"

    # Look for amino acid pattern
    # This SMARTS pattern looks for:
    # - A carbon with a carboxylic acid (-C(=O)OH)
    # - Connected to a nitrogen (can be NH2, NHR, or NR2)
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4]C(=O)[OH]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid group found"

    # Additional check for aromatic ring connected to amino acid portion
    # We'll look at the distance between aromatic atoms and the amino acid carbon
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if not amino_acid_matches:
        return False, "Could not locate amino acid group atoms"
        
    for match in amino_acid_matches:
        amino_c = match[1]  # The carbon connected to both NH2 and COOH
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic():
                # Check if aromatic ring is within reasonable distance of amino acid group
                path_length = len(Chem.GetShortestPath(mol, atom.GetIdx(), amino_c))
                if path_length <= 4:  # Within 3 bonds typically
                    return True, "Contains aromatic ring connected to amino acid group"

    return False, "Aromatic ring not properly connected to amino acid group"