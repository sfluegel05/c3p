"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    A glycosaminoglycan is characterized by polysaccharides containing substantial proportions of aminomonosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more generalized pattern for aminosugar (contains nitrogen in cyclic form)
    # Recognize aminosugar as any cyclic structure having nitrogen bonded to a carbon in the ring
    amino_cyclic_pattern = Chem.MolFromSmarts("C1[NH2,NH,N]C(O)C(O)C1")  # Generic pattern for aminocyclic structures
    if not mol.HasSubstructMatch(amino_cyclic_pattern):
        return False, "No aminomonosaccharide residues found"

    # Extending length consideration to qualify as substantial polysaccharide
    aminomonosaccharide_matches = mol.GetSubstructMatches(amino_cyclic_pattern)
    if len(aminomonosaccharide_matches) < 3:
        return False, "Does not meet the threshold for aminomonosaccharide residues"

    return True, "Contains polysaccharide chain with substantial aminomonosaccharide residues"