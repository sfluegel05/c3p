"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphoinositol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phospho-myo-inositol group:  Cyclohexane with five hydroxyls (1D-myo-inositol pattern)
    inositol_pattern = Chem.MolFromSmarts("C1(C(C(C(C(C1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring structure found"

    # Check for phosphate group attached to the inositol
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for glycerol backbone (C-C-C with at least one oxygen bonded)
    glycerol_pattern = Chem.MolFromSmarts("[C@H](CO)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Glycerophosphoinositols have long aliphatic chains often represented as carbon chains
    carbon_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCCC")
    chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if len(chain_matches) < 1:
        return False, "No adequate hydrophobic fatty acid chains detected"

    return True, "Contains essential glycerophosphoinositol structures with inositol ring, phosphate group, and attached fatty acid chains"

# Example usage:
# result, message = is_glycerophosphoinositol("P(O[C@H]1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCCCCCC/C=C\C/C=C\C/C=C\CCCCC)COC(=O)CCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC)(O)=O")
# print(result, message)