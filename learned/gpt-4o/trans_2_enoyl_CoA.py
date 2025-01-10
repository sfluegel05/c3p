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
    
    # Parse the SMILES into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the Coenzyme A core pattern
    coA_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@H]1OP(=O)(O)O)n2cnc3c(N)ncnc23")
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "Coenzyme A moiety is missing or incorrect"

    # Define pattern for the trans configuration of the alkene
    trans_double_bond_pattern = Chem.MolFromSmarts("C=C")  # Basic alkene pattern
    
    # Verify presence of trans double bond, using stereochemistry
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetStereo() in (Chem.rdchem.BondStereo.STEREOE, Chem.rdchem.BondStereo.STEREOTRANS):
            # Check that one end of this double bond leads to a thioester functionality
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            
            # Backtrack from alkene towards an thioester pattern: SC(=O)
            thioester_pattern = Chem.MolFromSmarts("SC(=O)")
            if atom1.HasSubstructMatch(thioester_pattern) or atom2.HasSubstructMatch(thioester_pattern):
                return True, "Detected Coenzyme A moiety with trans-2-enoyl moiety and appropriate thioester linkage"
    
    return False, "Trans-2-enoyl moiety with thioester linkage not found or improperly configured"