"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    A mucopolysaccharide is a polysaccharide composed of alternating units from uronic acids and glycosamines,
    and commonly partially esterified with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Redefine SMARTS pattern for uronic acid: carboxyl connected to carbohydrate
    uronic_acid_pattern = Chem.MolFromSmarts("C(=O)[O-]")

    # Redefine SMARTS pattern for glycosamine: presence of nitrogen in a sugar ring
    glycosamine_pattern = Chem.MolFromSmarts("[C,CX4]-[NX3][CX4]")

    # Redefine pattern for sulfuric acid ester: usually involves OS(=O)(=O) linkages to carbon
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[O,C]")

    # Check for presence of patterns
    uronic_acid_matches = mol.GetSubstructMatches(uronic_acid_pattern)
    glycosamine_matches = mol.GetSubstructMatches(glycosamine_pattern)
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)

    # If both uronic acid and glycosamine units are present
    if uronic_acid_matches and glycosamine_matches:
        # Check if there's adequate connectivity between these key units
        for ua in uronic_acid_matches:
            for ga in glycosamine_matches:
                # Find connectivity between matches
                if any(mol.GetBondBetweenAtoms(ua_atom, ga_atom) for ua_atom in ua for ga_atom in ga):
                    if sulfate_matches:  
                        return True, "Contains uronic acids and glycosamines with potential sulfate esterification"
                    else:
                        return True, "Contains uronic acids and glycosamines but no sulfuric acid esterification"

    return False, "No characteristic pattern of uronic acids and glycosamines found or insufficient connectivity"