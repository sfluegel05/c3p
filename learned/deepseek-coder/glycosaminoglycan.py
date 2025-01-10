"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    Glycosaminoglycans are polysaccharides containing a substantial proportion of aminomonosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible patterns for polysaccharide detection
    # Look for multiple sugar units connected by glycosidic bonds
    glycosidic_pattern = Chem.MolFromSmarts("[OX2;H0][CX4][CX4][OX2;H0]")  # O-C-C-O
    sugar_ring_pattern = Chem.MolFromSmarts("[CX4][CX4][OX2]")  # C-C-O in ring
    
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    sugar_ring_matches = mol.GetSubstructMatches(sugar_ring_pattern)
    
    if len(glycosidic_matches) < 1 or len(sugar_ring_matches) < 1:
        return False, "Insufficient glycosidic bonds or sugar rings for a polysaccharide"

    # More comprehensive amino sugar detection
    # Look for amino groups (-NH2, -NH-), N-acetyl groups, and sulfated amino groups
    amino_sugar_patterns = [
        Chem.MolFromSmarts("[CX4][NX3H2]"),  # -C-NH2
        Chem.MolFromSmarts("[CX4][NX3H][CX4]"),  # -C-NH-C-
        Chem.MolFromSmarts("[CX4][NX3]([CX3](=[OX1]))[CX4]"),  # N-acetyl
        Chem.MolFromSmarts("[CX4][NX3]([SX4](=[OX1])(=[OX1])[OX2-])[CX4]")  # Sulfated amino
    ]
    
    amino_sugar_count = 0
    for pattern in amino_sugar_patterns:
        amino_sugar_count += len(mol.GetSubstructMatches(pattern))
    
    if amino_sugar_count < 1:  # Require at least 1 amino sugar residue
        return False, "Insufficient amino sugar residues for a glycosaminoglycan"

    # Check molecular weight - glycosaminoglycans are typically large molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:  # Lowered threshold
        return False, "Molecular weight too low for a glycosaminoglycan"

    # Count carbons, oxygens, and nitrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 10:  # Lowered minimum
        return False, "Too few carbons for a glycosaminoglycan"
    if o_count < 5:  # Lowered minimum
        return False, "Too few oxygens for a glycosaminoglycan"
    if n_count < 1:  # Lowered minimum
        return False, "Too few nitrogens for a glycosaminoglycan"

    return True, "Contains polysaccharide structure with significant amino sugar residues"