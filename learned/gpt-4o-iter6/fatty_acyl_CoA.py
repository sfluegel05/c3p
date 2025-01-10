"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
from rdkit import Chem

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    A fatty acyl-CoA contains a Coenzyme A conjugated to a fatty acyl chain via a thioester linkage.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns should match with entire Coenzyme A structure, including phosphates and adenine.
    # Improved Coenzyme A pattern: adenosine structure + phosphates + thiol
    coenzyme_a_pattern = Chem.MolFromSmarts("O=C([C@@H](C(=O)NCCSC(=O)C))SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)N2C=NC3=C2N=CN=C3N")
    if not mol.HasSubstructMatch(coenzyme_a_pattern):
        return False, "Complete Coenzyme A structure not detected"
    
    # Thioester linkage pattern (R-C(=O)-S-CoA) focusing on correct acyl linkage
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)C")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage with CoA not found"
    
    # Fatty acyl chain should be detected based on hydrocarbon length and flexibility
    # Improved hydrocarbon chain pattern to allow more variation with potential double bonds
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("C(~C)(~[CH,CH2]){8,}")
    if not mol.HasSubstructMatch(hydrocarbon_chain_pattern):
        return False, "No sufficient hydrocarbon chain detected"

    return True, "Contains complete Coenzyme A structure, thioester linkage, and fatty acyl moiety"