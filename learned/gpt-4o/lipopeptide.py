"""
Classifies: CHEBI:46895 lipopeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide consists of a peptide with attached lipid moieties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Define SMARTS pattern for a peptide bond (-CO-NH-)
    peptide_bond_pattern = Chem.MolFromSmarts("[C](=O)[N]")
    
    if not mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "No peptide bonds found"

    # Check for lipid moieties: long hydrocarbon chains
    # Simple pattern to match long sequences of carbon atoms
    long_chain_pattern = Chem.MolFromSmarts("[CH2,CH3]~[CH2,CH3]~[CH2,CH3]~[CH2,CH3]~[CH2,CH3]~[CH2,CH3]~[CH2,CH3]~[CH2,CH3]")
    
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long hydrocarbon chains found"

    return True, "Contains both peptide bonds and long hydrocarbon chains indicative of lipids"

# Example usage:
# smiles_string = "CCCCC(=O)NCC(=O)NCCC(=O)OC1=C(NC2=C1C=CC=C2)O"
# is_lipopeptide(smiles_string)