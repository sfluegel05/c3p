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
        return (False, "Invalid SMILES string")

    # Define patterns for detecting glycosidic linkages
    # Adjust the SMARTS to capture possible glycosidic bonding configurations
    glycosidic_patterns = [
        Chem.MolFromSmarts("O[C@H]1[C@@H](O)O[C@@H](C)[C@H](CO)O1"),  # Generic pyranoside
        Chem.MolFromSmarts("C[C@H]1O[C@H](O[C@H]2O[C@H]([C@H](O)[C@@H]2O)CO)C(O)[C@H](O)[C@@H]1O"),  # Example potential linkage
        Chem.MolFromSmarts("C1(CO)C(O)C(O)O[C@@H]1[*]"),  # Cyclic ethers
        Chem.MolFromSmarts("C1(O)C(O[*])C(O)C(O)C(O)C1"),  # Another ether form
        Chem.MolFromSmarts("C([O])[C@H](O)[C@H](O)[C@H](O)CO"),  # Open form
    ]
    
    # Check for presence of glycosidic linkages
    has_glycosidic_linkage = any(mol.HasSubstructMatch(pattern) for pattern in glycosidic_patterns)
    if not has_glycosidic_linkage:
        return (False, "No glycosidic linkages found")
    
    # Better pattern for mono-saccharide residues like glucose
    saccharide_patterns = [
        Chem.MolFromSmarts("C1C(O)C(O)C(O)C(O)C1"),  # Simple glucosyl structure
        Chem.MolFromSmarts("[C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)O1"),  # More generalized
        Chem.MolFromSmarts("C1C(O)C(O)C2(O1)OCC(roval)[C@H](O)[C@H](O)C2"), # Furanosyl
    ]

    matches = set()
    for pattern in saccharide_patterns:
        matches.update(mol.GetSubstructMatches(pattern))
    
    if len(matches) <= 10:
        return (False, f"Only found {len(matches)} monosaccharide units, require more than 10")

    # Count carbon and oxygen
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Ensures high oxygen to carbon ratio
    if o_count / c_count < 0.5:
        return False, f"Low oxygen to carbon ratio: {o_count}/{c_count}"
    
    return (True, "Valid polysaccharide containing necessary glycosidic linkages and sufficient monosaccharide units")


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