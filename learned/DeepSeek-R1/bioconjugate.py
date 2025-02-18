"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate (contains at least two biological molecules linked covalently).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define SMARTS patterns for biological molecule parts
    amino_acid = Chem.MolFromSmarts("[NH2]-[CH]-[C](=O)[OH]")  # Basic amino acid pattern
    sugar = Chem.MolFromSmarts("[C]1O[C@H](O)[C@H](O)[C@H](O)[C@H]1O")  # Pyranose ring
    fatty_acid = Chem.MolFromSmarts("CCCCCCCCCCCC(=O)OH")  # Long chain carboxylic acid
    
    # Check for matches
    has_amino = mol.HasSubstructMatch(amino_acid)
    has_sugar = mol.HasSubstructMatch(sugar)
    has_fatty = mol.HasSubstructMatch(fatty_acid)
    
    # Count biological components
    bio_count = sum([has_amino, has_sugar, has_fatty])
    
    if bio_count >= 2:
        return True, "Contains multiple biological components"
    
    # Check for peptide bonds (amide links)
    amide = Chem.MolFromSmarts("[CX3](=O)[NX3H]")
    if len(mol.GetSubstructMatches(amide)) >= 2:
        return True, "Contains peptide bonds"
    
    # Check for disulfide bonds (S-S)
    disulfide = Chem.MolFromSmarts("[S][S]")
    if mol.HasSubstructMatch(disulfide):
        return True, "Contains disulfide bond"
    
    # Check for thioether (S connected to two carbons)
    thioether = Chem.MolFromSmarts("[S]([#6])[#6]")
    if mol.HasSubstructMatch(thioether):
        # Check if connected to a peptide or amino acid
        if has_amino or len(mol.GetSubstructMatches(amide)) >= 1:
            return True, "Thioether linked to peptide/amino acid"
    
    return False, "Does not meet bioconjugate criteria"