"""
Classifies: CHEBI:46895 lipopeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide consists of a peptide with an attached lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for peptide presence (N-C(=O)-C pattern in backbone)
    peptide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[C]")
    if not mol.HasSubstructMatch(peptide_pattern):
        return False, "No peptide backbone found"

    # Check for lipid presence more comprehensively (long hydrocarbon chains, possibly branched or cyclic)
    # Lipid pattern: Include flexible SMARTS to identify both linear and branched large carbons chains
    lipid_pattern = Chem.MolFromSmarts("[R0;!$(*=*)&!D1]-[R0;!$(*=*)&!D1]-[R0;!$(*=*)&!D1]-[R0;!$(*=*)&!D1]-[R0;!$(*=*)&!D1]-[R0;!$(*=*)&!D1]")
    
    for i in range(len(mol.GetSubstructMatches(lipid_pattern))):
        path = mol.GetSubstructMatch(lipid_pattern)
        if len(path) >= 8:
            return True, "Contains both peptide and lipid components"

    return False, "No recognizably complex lipid part found"