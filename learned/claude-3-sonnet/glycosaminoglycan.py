"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
"""
Classifies: CHEBI:36976 glycosaminoglycan
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    A glycosaminoglycan is a polysaccharide containing a substantial proportion of aminomonosaccharide residues.

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
    
    # Look for aminomonosaccharide residues
    aminosugar_pattern = Chem.MolFromSmarts("[NX3H2][CR]([OR])([OR])C(=O)[OR]")
    aminosugar_matches = mol.GetSubstructMatches(aminosugar_pattern)
    if len(aminosugar_matches) == 0:
        return False, "No aminomonosaccharide residues found"
    
    # Look for polysaccharide backbone
    backbone_pattern = Chem.MolFromSmarts("[OX2]C([OX2])[OX2]C([OX2])[OX2]C([OX2])[OX2]")
    backbone_matches = mol.GetSubstructMatches(backbone_pattern)
    if len(backbone_matches) == 0:
        return False, "No polysaccharide backbone found"
    
    # Check for characteristic glycosidic linkages
    glycosidic_pattern = Chem.MolFromSmarts("[OX2]C[OX2]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if len(glycosidic_matches) < 3:
        return False, "Insufficient glycosidic linkages for a glycosaminoglycan"
    
    # Calculate molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    
    # Check if properties are consistent with glycosaminoglycans
    if mol_wt < 500 or n_rotatable < 10 or n_rings < 5:
        return False, "Molecular properties not consistent with glycosaminoglycans"
    
    return True, "Contains aminomonosaccharide residues connected to a polysaccharide backbone"