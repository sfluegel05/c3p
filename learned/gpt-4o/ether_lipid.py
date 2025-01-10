"""
Classifies: CHEBI:64611 ether lipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid has one or more carbon atoms on the glycerol backbone linked
    to an alkyl chain via an ether linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a glycerol backbone pattern - 3 carbons in a row
    glycerol_pattern = Chem.MolFromSmarts("[CH2][CH][CH2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for ether linkage patterns - C-O-C
    ether_pattern = Chem.MolFromSmarts("C-O-C")
    if not mol.HasSubstructMatch(ether_pattern):
        return False, "No ether linkage found"

    # Ensure the presence of ether linkages over ester
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ether_matches) <= len(ester_matches):
        return False, "Not enough ether linkages compared to ester linkages"

    # Optional: Check for typical lipid components like phosphocholine
    phospho_pattern = Chem.MolFromSmarts("P(=O)([O-])(O)")
    if not mol.HasSubstructMatch(phospho_pattern):
        return False, "Missing phosphate group, possibly not a lipid"

    return True, "Ether lipid identified by ether linkages on glycerol backbone and presence of phosphate group."