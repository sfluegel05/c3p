"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
"""
Classifies: CHEBI:18086 glycosaminoglycan
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    A glycosaminoglycan is a polysaccharide containing a substantial proportion
    of aminomonosaccharide residues.

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
    
    # Check for aminomonosaccharide residues
    aminomonosaccharide_pattern = Chem.MolFromSmarts("[NX3+0;R]")
    aminomonosaccharide_matches = mol.GetSubstructMatches(aminomonosaccharide_pattern)
    if len(aminomonosaccharide_matches) < 2:
        return False, "Insufficient aminomonosaccharide residues"
    
    # Check for polysaccharide backbone
    polysaccharide_pattern = Chem.MolFromSmarts("[OX2]~[CX4]~[OX2]~[CX4]~[OX2]~[CX4]")
    polysaccharide_matches = mol.GetSubstructMatches(polysaccharide_pattern)
    if len(polysaccharide_matches) < 1:
        return False, "No polysaccharide backbone found"
    
    # Check for glycosidic bonds
    glycosidic_pattern = Chem.MolFromSmarts("[OX2]~[CX4]~[OX2]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if len(glycosidic_matches) < 2:
        return False, "Insufficient glycosidic bonds"
    
    # Check for sulfate or carboxylate groups (common in glycosaminoglycans)
    sulfate_pattern = Chem.MolFromSmarts("[O-;X1-]=[S+](=O)[O-]")
    carboxylate_pattern = Chem.MolFromSmarts("[O-;X1-]C(=O)[O-]")
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(sulfate_matches) + len(carboxylate_matches) < 1:
        return False, "No sulfate or carboxylate groups found"
    
    # Check molecular weight (glycosaminoglycans typically >500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for glycosaminoglycan"
    
    return True, "Contains aminomonosaccharide residues, polysaccharide backbone, glycosidic bonds, and sulfate/carboxylate groups"