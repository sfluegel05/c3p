"""
Classifies: CHEBI:28892 ganglioside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    A ganglioside is a glycosphingolipid containing one or more sialic acids 
    linked to an oligosaccharide chain attached to a ceramide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight - gangliosides are large molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800:
        return False, "Molecular weight too low for ganglioside"

    # Look for ceramide backbone (sphingosine + fatty acid)
    # Sphingosine base with long chain and OH/NH2 groups
    sphingosine_pattern = Chem.MolFromSmarts("[CH2]([OH,O])[CH]([NH,N])[CH]([OH,O])/C=C/CCCCCC")
    
    # Fatty acid chain attached to NH
    fatty_acid_pattern = Chem.MolFromSmarts("[NH]C(=O)CCCCCCC")
    
    if not (mol.HasSubstructMatch(sphingosine_pattern) and mol.HasSubstructMatch(fatty_acid_pattern)):
        return False, "No ceramide backbone found"

    # More precise sialic acid (Neu5Ac) pattern
    # Includes pyranose ring with carboxyl, N-acetyl, and specific OH positions
    sialic_pattern = Chem.MolFromSmarts("[C,c]([OH,O])[C,c]1[C,c]([C,c]([C,c]([C,c]([C,c]1[C,c](=O)[O;H,-])[OH,O])[NH,N]C(=O)C)[OH,O])[OH,O]")
    
    if sialic_pattern is None:
        return False, "Error in sialic acid pattern"
        
    sialic_matches = len(mol.GetSubstructMatches(sialic_pattern))
    if sialic_matches < 1:
        return False, "No sialic acid (Neu5Ac) found"

    # Look for sugar units (pyranose rings)
    # More specific pattern for common sugars in gangliosides (Glc, Gal, GalNAc)
    sugar_pattern = Chem.MolFromSmarts("[C,c]1[O,o][C,c]([C,c]([C,c]([C,c]1)[OH,O])[OH,O])[CH2][OH,O]")
    
    if sugar_pattern is None:
        return False, "Error in sugar pattern"
        
    sugar_matches = len(mol.GetSubstructMatches(sugar_pattern))
    if sugar_matches < 3:
        return False, "Insufficient oligosaccharide chain"

    # Check for glycosidic linkages between sugars
    glycosidic_pattern = Chem.MolFromSmarts("[C,c]1[O,o][C,c][C,c][C,c][C,c]1[O,o][C,c]")
    if glycosidic_pattern is None:
        return False, "Error in glycosidic pattern"
        
    glycosidic_matches = len(mol.GetSubstructMatches(glycosidic_pattern))
    if glycosidic_matches < 2:
        return False, "Insufficient glycosidic linkages"

    # Count key atoms to verify overall composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 30:
        return False, "Too few carbons for ganglioside structure"
    if o_count < 15:
        return False, "Too few oxygens for ganglioside structure"
    if n_count < 2:
        return False, "Too few nitrogens for ganglioside structure"

    return True, f"Contains ceramide backbone, {sialic_matches} sialic acid(s), and oligosaccharide chain with {glycosidic_matches} glycosidic linkages"