"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem

def is_polysaccharide(smiles):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide consists of more than ten monosaccharide residues linked glycosidically.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Define a more comprehensive pattern for monosaccharide units
    # Revised patterns considering aldohexoses and ketohexoses general structures
    hexose_pattern = Chem.MolFromSmarts("C1(O)[C&H=O&X4][C&H&X2][C&H&X2][C&H&X2][C&H&X2]O1") 
    derivative_pattern = Chem.MolFromSmarts("NC(=O)C")  

    # Count hexose units, but incorporate variability detection (typical polysaccharide chain residues)
    hexose_count = len(mol.GetSubstructMatches(hexose_pattern))
    derivative_count = len(mol.GetSubstructMatches(derivative_pattern))

    # Assess total identifiable sugar-like units
    total_count = hexose_count + derivative_count

    if total_count < 10:
        return False, f"Detected {total_count} monosaccharide-like residues, which is below the polysaccharide threshold."

    # Define pattern for glycosidic linkages
    glyco_linkage_pattern = Chem.MolFromSmarts("[C&H]O[C&H]")  # Simplified ether linkage
    glyco_link_count = len(mol.GetSubstructMatches(glyco_linkage_pattern))
    
    # Verify there are adequate bonds to form a polysaccharide chain
    if glyco_link_count < total_count - 1:
        return False, f"Insufficient glycosidic linkages ({glyco_link_count}) compared to expected ({total_count - 1})"

    # Passed all checks
    return True, "Sufficient monosaccharide units and linkages to qualify as a polysaccharide"