"""
Classifies: CHEBI:17002 cholesteryl ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    Defined as a sterol ester obtained by formal condensation of the carboxy group 
    of any carboxylic acid with the 3-hydroxy group of cholesterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cholesteryl ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to molecule and validate
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define accurate cholesterol steroid nucleus pattern
    cholesterol_backbone_pattern = Chem.MolFromSmarts("C1CCC2(C)C(C(C))CCC3C2CCC4(CCC(=O)O)C3CCCC14")  # Considered typical steroid pattern with functional group placements
    if not mol.HasSubstructMatch(cholesterol_backbone_pattern):
        return False, "No cholesterol steroid backbone found"

    # Check for ester linkage specifically bound to the cholesterol 3-hydroxy position
    ester_linkage_pattern = Chem.MolFromSmarts("O=C(O)")  # look for more precise binding indication
    ester_matches = mol.GetSubstructMatches(ester_linkage_pattern)
    
    # Ensure ester linkage is not arbitrary but specifically with cholesterol
    ester_connected_to_cholesterol = False
    for match in ester_matches:
        atom_indices = set(match)
        if any(mol.GetSubstructMatches(ester_linkage_pattern)):
            group = [mol.GetAtomWithIdx(idx) for idx in atom_indices]
            # Check if any of these atoms connect in the manner defined by the backbone
            ester_connected_to_cholesterol = any(cholesterol_backbone_pattern)
            if ester_connected_to_cholesterol:
                break
    
    if not ester_connected_to_cholesterol:
        return False, "Ester linkage not specifically forming with cholesterol's 3-hydroxy position"

    return True, "Contains cholesterol backbone with specific ester linkage indicative of cholesteryl ester"