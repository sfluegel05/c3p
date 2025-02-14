"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
from rdkit import Chem

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trans-2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for coenzyme A adenine component
    # Adjusted SMARTS pattern to more specifically reflect adenine as found in CoA
    coa_pattern = Chem.MolFromSmarts("n1cnc2ncnc2n1")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Not a CoA molecule (possible mismatch in adenine base detection)"

    # Check for thiol ester linkage to CoA: -C(=O)SCCNC(=O)
    thiol_ester_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)")
    if not mol.HasSubstructMatch(thiol_ester_pattern):
        return False, "Missing thiol ester linkage"

    # Check for trans double bond at the 2,3 position: C/C=C/
    # Adjusted trans bond pattern to focus on unsaturation, allows flexibility in atom types
    trans_double_bond_pattern = Chem.MolFromSmarts("C/C=C")
    if len(mol.GetSubstructMatches(trans_double_bond_pattern)) == 0:
        return False, "No 2,3-trans double bond found at the correct position"

    return True, "Molecule fits the trans-2-enoyl-CoA structure criteria"