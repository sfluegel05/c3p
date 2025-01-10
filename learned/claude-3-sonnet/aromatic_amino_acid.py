"""
Classifies: CHEBI:33856 aromatic amino acid
"""
"""
Classifies: aromatic amino acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdDecomposition import IsAromatic

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
    for ring in mol.GetRingInfo().AtomRings():
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            has_aromatic = True
            break
            
    if not has_aromatic:
        return False, "No aromatic ring found"

    # Look for carboxylic acid group (-COOH)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"

    # Look for amino group (-NH2) or substituted amino (-NHR, -NR2)
    # We'll check for both primary amines and substituted amines
    primary_amine_pattern = Chem.MolFromSmarts("[NX3H2]")
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3H1]")
    tertiary_amine_pattern = Chem.MolFromSmarts("[NX3H0]")
    
    has_amine = (mol.HasSubstructMatch(primary_amine_pattern) or 
                 mol.HasSubstructMatch(secondary_amine_pattern) or
                 mol.HasSubstructMatch(tertiary_amine_pattern))
    
    if not has_amine:
        return False, "No amino group found"

    # Additional check to ensure the molecule has the basic characteristics
    # Count carbons, nitrogens and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 3:
        return False, "Too few carbons for an aromatic amino acid"
    if n_count < 1:
        return False, "Must contain at least one nitrogen"
    if o_count < 2:
        return False, "Must contain at least two oxygens (carboxylic acid)"

    return True, "Contains both aromatic ring and amino acid functionality"