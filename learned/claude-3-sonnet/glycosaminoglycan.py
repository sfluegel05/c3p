"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
"""
Classifies: glycosaminoglycan
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    
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
        
    # Count atoms and check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:  # GAGs are typically large molecules
        return False, "Molecular weight too low for glycosaminoglycan"

    # Count key atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if o_count < 5:  # Need multiple oxygen atoms for sugar rings
        return False, "Too few oxygen atoms for polysaccharide structure"
    
    if n_count == 0:  # Need nitrogen for amino sugars
        return False, "No nitrogen atoms found - required for aminosugar"

    # Look for pyranose/furanose ring patterns (sugar rings)
    sugar_pattern = Chem.MolFromSmarts("[CR1]1[CR1][CR1][CR1][CR1]O1")  # pyranose ring
    sugar_matches = len(mol.GetSubstructMatches(sugar_pattern))
    
    if sugar_matches == 0:
        return False, "No sugar rings detected"
        
    # Look for amino sugar pattern (-NH-CH-)
    amino_sugar = Chem.MolFromSmarts("[NX3H][CH1,CH2][OH1,OR]")
    amino_sugar_matches = len(mol.GetSubstructMatches(amino_sugar))
    
    if amino_sugar_matches == 0:
        return False, "No amino sugar residues detected"

    # Look for glycosidic linkages (-O- between sugars)
    glycosidic = Chem.MolFromSmarts("[CR1]O[CR1]")
    glycosidic_matches = len(mol.GetSubstructMatches(glycosidic))
    
    if glycosidic_matches < 1:
        return False, "No glycosidic linkages found"

    # Optional: Check for common modifications
    sulfate = Chem.MolFromSmarts("OS(=O)(=O)[OH1,O-]")
    carboxyl = Chem.MolFromSmarts("C(=O)[OH1,O-]")
    has_sulfate = mol.HasSubstructMatch(sulfate)
    has_carboxyl = mol.HasSubstructMatch(carboxyl)

    # Construct reason string
    features = []
    features.append(f"Contains {sugar_matches} sugar rings")
    features.append(f"{amino_sugar_matches} amino sugar residues")
    features.append(f"{glycosidic_matches} glycosidic linkages")
    if has_sulfate:
        features.append("sulfate groups")
    if has_carboxyl:
        features.append("carboxyl groups")
    
    reason = "Classified as glycosaminoglycan: " + ", ".join(features)

    # Final classification
    # Must have multiple sugar rings, amino sugars, and glycosidic linkages
    is_gag = (sugar_matches >= 1 and 
              amino_sugar_matches >= 1 and 
              glycosidic_matches >= 1)

    return is_gag, reason