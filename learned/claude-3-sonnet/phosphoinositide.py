"""
Classifies: CHEBI:18179 phosphoinositide
"""
"""
Classifies: CHEBI:18035 phosphoinositide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol that is phosphorylated at one or more of the
    hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for inositol ring attached to a phosphate group
    inositol_pattern = Chem.MolFromSmarts("OC1C(O)C(O)C(O)C(O)C1OP")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring attached to a phosphate group"
    
    # Check for additional phosphate groups on inositol
    inositol_p_pattern = Chem.MolFromSmarts("OC1C(O)C(OP)C(O)C(O)C1")
    if not mol.HasSubstructMatch(inositol_p_pattern):
        return False, "Inositol ring not phosphorylated"

    # Check molecular weight - phosphoinositides typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for phosphoinositide"

    # Count phosphorus atoms - phosphoinositides should have at least one
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count < 1:
        return False, "No phosphorus atoms found"

    return True, "Contains an inositol ring attached to a phosphate group, with additional phosphate groups on the inositol ring"