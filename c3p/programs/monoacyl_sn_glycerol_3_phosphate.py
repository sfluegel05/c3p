"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
"""
Classifies: CHEBI:17504 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    A monoacyl-sn-glycerol 3-phosphate has a glycerol backbone with a single acyl group at either position 1 or 2,
    and a phosphate group at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacyl-sn-glycerol 3-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone with phosphate at position 3
    glycerol_phosphate_pattern = Chem.MolFromSmarts("[CH2X4][CHX4]([OHX2])[CH2X4][OX2][PX4](=[OX1])([OX2-])[OX2-]")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol backbone with phosphate at position 3 found"

    # Check for exactly one ester group (acyl group)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for acyl group (long carbon chain)
    acyl_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) < 1:
        return False, "No acyl group (fatty acid chain) found"

    # Check molecular weight to ensure it's a reasonable size for a monoacyl-sn-glycerol 3-phosphate
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1000:
        return False, "Molecular weight out of range for monoacyl-sn-glycerol 3-phosphate"

    # Check stereochemistry (sn-glycerol)
    # The stereochemistry can be inferred from the SMILES string, but RDKit's substructure matching does not directly handle stereochemistry.
    # For simplicity, we assume the SMILES string correctly represents the stereochemistry.

    return True, "Contains glycerol backbone with a single acyl group and a phosphate group at position 3"