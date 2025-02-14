"""
Classifies: CHEBI:36500 glucosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    A glucosylceramide contains a sphingosine or sphinganine backbone, a fatty acid attached via an amide bond and a beta-D-glucose head group linked via the 1 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosylceramide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for terminal beta-D-glucose linked via the 1 position
    glucose_pattern = Chem.MolFromSmarts("[OX2][C@H]1[C@@H]([C@@H](O)[C@H](O)[C@H](CO)O)O[C@H]1[CX4]")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No terminal beta-D-glucose head group linked at the 1 position found"

    # Check for sphingosine/sphinganine backbone and amide linkage
    sphingosine_pattern = Chem.MolFromSmarts("[C@H]([OX2])[C@H]([NX3][CX3](=[OX1]))[C@H]([OX2])")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine or sphinganine backbone with amide linkage found"
    
    #Check for amide bond to the sphingosine nitrogen
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3][C@H]([OX2])")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found linked to the sphingosine backbone."
    
    # Check for long chain fatty acid (at least 8 carbons, can be more sophisticated)
    
    amide_match = mol.GetSubstructMatches(amide_pattern)
    if amide_match:
        for match in amide_match:
            amide_carbon = match[0]
            fatty_acid_carbons = 0
            for atom in mol.GetAtoms():
                if atom.GetIdx() == amide_carbon:
                    for neighbor in atom.GetNeighbors():
                        fatty_acid_carbons = sum(1 for a in Chem.GetShortestPath(mol,neighbor.GetIdx(),amide_carbon) if a.GetAtomicNum() == 6)
                        if fatty_acid_carbons < 7:
                            return False, "No long chain fatty acid found"

    return True, "Contains a beta-D-glucose head group, sphingosine/sphinganine backbone and a fatty acid via an amide bond"