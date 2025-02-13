"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate includes a polyprenol chain attached via a phosphate ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Recognize phosphate (P=O with single or two oxygens linked) or diphosphate groups
    phosphate_pattern = Chem.MolFromSmarts("O=P([O-])([O-])" )  # Adjusted for phosphorylated oxygen
    debug_info = []

    oid = mol.HasSubstructMatch(phosphate_pattern)
    debug_info.append(f"Phosphate group match: {oid}")

    if not oid:
        return False, "No phosphate or diphosphate group found or incorrectly connected"
    
    # Improved ester linkage search targeting polyprenol specific attachment
    ester_linkage = Chem.MolFromSmarts("O=P([O-])([O-])O")
    match_ester = mol.HasSubstructMatch(ester_linkage)
    debug_info.append(f"Ester linkage to phosphate match: {match_ester}")

    if not match_ester:
        return False, "Ester linkage to phosphate group not found"

    # Refine isoprene detection, ensuring correct context inside polyprenyl chain.
    isoprene_pat = Chem.MolFromSmarts("C(C)=C(C)CC")  # Aligning more closely to polyprenol structure implications
    isoprene_match = len(mol.GetSubstructMatches(isoprene_pat))

    debug_info.append(f"Isoprene unit counts: {isoprene_match}")

    if isoprene_match < 3:  # Typically three or more isoprene units
        return False, f"Insufficient isoprene units for polyprenol chain: {isoprene_match}"
    
    # Optionally add other layers of validation here if needed.

    return True, "Successfully identified as a polyprenol phosphate" + "; Debug info: " + ", ".join(debug_info)