"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: CHEBI:34555 chalcone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_chalcone(smiles: str):
    """
    Determines if a molecule is a chalcone or chalcone derivative based on its SMILES string.
    A chalcone is a ketone with the structure ArCH=CH(=O)Ar, where Ar represents an aryl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for chalcone backbone pattern (ArCH=CH(=O)Ar)
    chalcone_pattern = Chem.MolFromSmarts("[c:1]=[c:2][CH:3]=[CH:4][C:5]=[O:6].[c:7]")
    matches = mol.GetSubstructMatches(chalcone_pattern)
    if not matches:
        return False, "Chalcone backbone pattern not found"

    # Check for aromaticity in the rings
    rings = mol.GetRingInfo().AtomRings()
    aromatic_rings = [ring for ring in rings if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]
    if len(aromatic_rings) < 2:
        return False, "At least two aromatic rings are required for chalcones"

    # Check for common chalcone substituents
    substituent_patterns = [
        Chem.MolFromSmarts("[OH]"),  # hydroxyl
        Chem.MolFromSmarts("[OC]"),  # methoxy
        Chem.MolFromSmarts("[CC=C(C)C]"),  # prenyl
        Chem.MolFromSmarts("[F,Cl,Br,I]"),  # halides
        Chem.MolFromSmarts("[OC(=O)C]"),  # ester
        Chem.MolFromSmarts("[C(=O)N]"),  # amide
        Chem.MolFromSmarts("[C(=O)O]"),  # carboxyl
        Chem.MolFromSmarts("[N+]"),  # quaternary amines
        Chem.MolFromSmarts("[O-]"),  # phenoxide
    ]
    has_substituents = any(mol.HasSubstructMatch(pattern) for pattern in substituent_patterns)

    # Return the classification result
    if has_substituents:
        return True, "Contains chalcone backbone and common substituents"
    else:
        return True, "Contains chalcone backbone"