"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    A glycerophosphoinositol is a glycerophospholipid with inositol esterified to the phosphate group at sn-3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool, str: True with reason if molecule is a glycerophosphoinositol, False otherwise with reason
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycerol backbone pattern: C-C-C with oxygens
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)COP(=O)(O)O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with phosphate at sn-3 found"
    
    # Inositol moiety attached to phosphate
    inositol_phosphate_pattern = Chem.MolFromSmarts("P(=O)(OCC1OC(O)C(O)C(O)C1O)(OC)O")
    if not mol.HasSubstructMatch(inositol_phosphate_pattern):
        return False, "No inositol group esterified to phosphate found"
    
    # Fatty acyl chains pattern: detect at least two ester linkages typical for glycerophospholipids
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Expected at least 2 ester linkages, found {len(ester_matches)}"

    return True, "Contains glycerol backbone with inositol esterified to phosphate at sn-3 position"

# Example usage: 
# smiles_input = "P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCC)(O)=O"
# is_glycerophosphoinositol(smiles_input)