"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
"""
Classifies: polyprenol phosphate
Definition: A prenol phosphate resulting from the formal condensation of the terminal allylic hydroxy group 
of a polyprenol with 1 mol eq. of phosphoric acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Look for phosphate group(s)
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])[OX2]")
    diphosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2])OP(=O)([OX2])[OX2]")
    
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)
    has_diphosphate = mol.HasSubstructMatch(diphosphate_pattern)
    
    if not (has_phosphate or has_diphosphate):
        return False, "No phosphate/diphosphate group found"

    # Look for isoprene units (C=C-C(C)-C)
    isoprene_pattern = Chem.MolFromSmarts("C=C-C(C)-C")
    isoprene_matches = len(mol.GetSubstructMatches(isoprene_pattern))
    
    if isoprene_matches < 2:
        return False, "Not enough isoprene units for a polyprenol"
    
    # Check for allylic oxygen-phosphorus connection
    allylic_phosphate = Chem.MolFromSmarts("C=C-COP(=O)")
    if not mol.HasSubstructMatch(allylic_phosphate):
        return False, "Phosphate not connected to allylic position"
        
    # Count carbons to verify chain length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:  # Minimum 2 isoprene units (C10)
        return False, "Carbon chain too short for polyprenol"
        
    # Look for characteristic branching pattern of polyprenols
    methyl_branches = Chem.MolFromSmarts("C-C(C)=C")
    branch_count = len(mol.GetSubstructMatches(methyl_branches))
    
    if branch_count < 2:
        return False, "Not enough methyl branches for polyprenol structure"
        
    # Check for conjugated double bonds pattern
    conjugated_pattern = Chem.MolFromSmarts("C=CC=C")
    if mol.HasSubstructMatch(conjugated_pattern):
        return False, "Contains conjugated double bonds, not characteristic of polyprenols"
        
    # Calculate number of rotatable bonds to verify flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Not enough rotatable bonds for polyprenol chain"

    # Success case - provide detailed reason
    phosphate_type = "diphosphate" if has_diphosphate else "phosphate"
    return True, f"Contains polyprenol chain ({isoprene_matches} isoprene units) with {phosphate_type} group at allylic position"