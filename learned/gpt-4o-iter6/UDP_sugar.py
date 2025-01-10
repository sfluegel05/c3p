"""
Classifies: CHEBI:17297 UDP-sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    A UDP-sugar is a pyrimidine nucleotide-sugar having UDP as the nucleotide component attached
    to an unspecified sugar via an anomeric diphosphate linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an UDP-sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for Uridine (nucleoside part of UDP)
    uridine_pattern = Chem.MolFromSmarts("n1ccc(=O)[nH]c1=O")
    if not mol.HasSubstructMatch(uridine_pattern):
        return False, "No uridine found"

    # SMARTS pattern for diphosphate link
    diphosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)OP(O)(=O)")
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "No diphosphate linkage found"

    # Look for sugar moiety attached via diphosphate (unrestricted sugar form)
    sugar_attachment_pattern = Chem.MolFromSmarts("O[C@H]1[C@H]([C@@H]([C@H]([C@H]1O)O)O)OP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2ccc(=O)[nH]c2=O")
    if not mol.HasSubstructMatch(sugar_attachment_pattern):
        return False, "No suitable sugar attachment found via diphosphate linkage"

    return True, "Contains UDP group with appropriate sugar attachment"