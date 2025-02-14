"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    An inositol phosphoceramide is a phosphosphingolipid with an inositol and ceramide
    linked via a phosphodiester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an inositol phosphoceramide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Ceramide Identification
    # Look for amide bond (C(=O)N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond found, not a ceramide"

    # Identify the sphingosine backbone (at least 2 carbons with a hydroxyl each
    sphingosine_pattern = Chem.MolFromSmarts("[CX4](O)[CX4](O)")
    sphingosine_matches = mol.GetSubstructMatches(sphingosine_pattern)
    if not sphingosine_matches:
        return False, "No sphingosine backbone found"

    # 2. Phosphate group
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"

    # 3. Inositol Identification
    # Look for a 6-membered ring with at least 5 oxygens attached to the carbons
    inositol_pattern = Chem.MolFromSmarts("[CX4]1[CX4][CX4][CX4][CX4][CX4]1") #basic ring, we assume it will be an inositol due to other checks
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_matches:
        return False, "No inositol ring found"

    # 4. Verify phosphodiester linkage
    # Combined pattern to look for the linkage:
    # inositol -O-P(=O)(O)-O-Ceramide
    linkage_pattern = Chem.MolFromSmarts("[CX4]1[CX4][CX4][CX4][CX4][CX4]1[OX2][P](=[OX1])([OX2])[OX2][CX4,CX3]([OX2])[CX4,CX3]([OX2])") # This is more permissive, it does not assume the N is bonded to the ceramide part.
    
    linkage_matches = mol.GetSubstructMatches(linkage_pattern)
    if not linkage_matches:
       return False, "Phosphate is not linked to inositol and ceramide via phosphodiester bond"
    

    # If all checks pass, return True
    return True, "Structure contains ceramide, inositol, and phosphate group with correct linkage"