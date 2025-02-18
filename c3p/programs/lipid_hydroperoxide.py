"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is characterized by the presence of hydroperoxy groups (-OOH) 
    attached to a lipid-like structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for hydroperoxy group pattern (-OOH)
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2H]O")
    hydroperoxy_matches = mol.GetSubstructMatches(hydroperoxy_pattern)
    if not hydroperoxy_matches:
        return False, "No hydroperoxy groups found"
    
    # Count carbon atoms to verify lipid-like structure (heuristic: must have at least 12 carbon atoms)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12:
        return False, f"Only {c_count} carbon atoms found, not enough for typical lipid structure"

    # Check for presence of carboxylic group or long-chain fatty-acid like structures
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_pattern)

    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not carboxylic_matches and not ester_matches:
        return False, "No carboxylic or long-chain ester structures found"

    return True, f"Contains {len(hydroperoxy_matches)} hydroperoxy group(s) and lipid-like structure with {c_count} carbon atoms"