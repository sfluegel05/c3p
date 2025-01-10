"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide is composed of more than ten monosaccharide residues linked glycosidically.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycosidic linkages
    # This is an overly simplified pattern, but it captures the essence
    glycosidic_pattern = Chem.MolFromSmarts("[O&r3]([C&r3])[C&r3]")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkages found"
    
    # Count potential mono- and di-saccharide units (e.g., glucose-like)
    monosaccharide_pattern = Chem.MolFromSmarts("C1C(O)C(O)C(O)C(O)C1")
    matches = mol.GetSubstructMatches(monosaccharide_pattern)
    
    if len(matches) <= 10:
        return False, f"Only found {len(matches)} monosaccharide units, require more than 10"

    # Count carbohydrate-specific features
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Typical polysaccharides have high O/C ratios
    if o_count / c_count < 0.5:
        return False, f"Low oxygen to carbon ratio: {o_count}/{c_count}"

    return True, "Contains glycosidic linkages and more than 10 monosaccharide-like units"

__metadata__ = { 
    'chemical_class': { 
        'id': 'CHEBI:16646',
        'name': 'polysaccharide',
        'definition': 'A biomacromolecule consisting of large numbers of monosaccharide residues linked glycosidically. Includes those containing more than ten monosaccharide residues.',
        'parents': ['CHEBI:16895', 'CHEBI:18154']
    },
    'config': {
        'test_proportion': 0.1
    }
}