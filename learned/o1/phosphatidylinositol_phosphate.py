"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
"""
Classifies: phosphatidylinositol phosphate
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    A phosphatidylinositol phosphate has a glycerol backbone with two fatty acid chains,
    a phosphate group at position 3 linked to an inositol ring, and additional phosphate
    groups on the inositol ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylinositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Suppress sanitization warnings
    Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE)
    
    # Identify glycerol backbone with two ester-linked fatty acids
    glycerol_fatty_acid_pattern = Chem.MolFromSmarts("""
    [C@H]([O][C](=O)[C,C])  # sn-1 ester
    [C@@H]([O][C](=O)[C,C]) # sn-2 ester
    [C@H]([O][P](=O)([O])[O]) # sn-3 phosphate
    """)
    if not mol.HasSubstructMatch(glycerol_fatty_acid_pattern):
        return False, "Glycerol backbone with two fatty acid chains not found"
    
    # Identify inositol ring attached via phosphodiester bond
    inositol_phosphate_pattern = Chem.MolFromSmarts("""
    [O][P](=O)([O])[O][C@H]1[C@H]([O])[C@H]([O])[C@@H]([O])[C@H]([O])[C@@H]1[O] # myo-inositol ring
    """)
    if not mol.HasSubstructMatch(inositol_phosphate_pattern):
        return False, "Inositol ring attached via phosphodiester linkage not found"
    
    # Identify additional phosphate groups on inositol ring
    # Phosphate groups attached to inositol hydroxyls
    inositol_additional_phosphate_pattern = Chem.MolFromSmarts("""
    [C@H]1([O][P](=O)([O])[O])[C@H]([O])[C@@H]([O])[C@H]([O])[C@@H]([O])[C@H]1[O] |
    [C@H]1([O])[C@H]([O][P](=O)([O])[O])[C@@H]([O])[C@H]([O])[C@@H]([O])[C@H]1[O] |
    [C@H]1([O])[C@H]([O])[C@@H]([O][P](=O)([O])[O])[C@H]([O])[C@@H]([O])[C@H]1[O] |
    [C@H]1([O])[C@H]([O])[C@@H]([O])[C@H]([O][P](=O)([O])[O])[C@@H]([O])[C@H]1[O] |
    [C@H]1([O])[C@H]([O])[C@@H]([O])[C@H]([O])[C@@H]([O][P](=O)([O])[O])[C@H]1[O]
    """)
    if not mol.HasSubstructMatch(inositol_additional_phosphate_pattern):
        return False, "Additional phosphate groups on inositol ring not found"
    
    return True, "Molecule is a phosphatidylinositol phosphate"