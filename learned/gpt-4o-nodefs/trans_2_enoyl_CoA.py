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

    # Look for trans double bond - more generic pattern targeting common trans enoyl structure
    trans_double_bond_patterns = [
        Chem.MolFromSmarts("C/C=C/C"),
        Chem.MolFromSmarts("C\\C=C\\C")
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in trans_double_bond_patterns):
        return False, "No trans double bond found"
    
    # Look for key parts of the CoA structure - divided into smaller components
    coa_substructures = [
        Chem.MolFromSmarts("C(=O)SCCNC(=O)"),  # Part of the CoA linkage
        Chem.MolFromSmarts("C(=O)[C@H](O)C(C)(C)COP(O)(O)=O"),  # Phosphopantetheine
        Chem.MolFromSmarts("OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1")  # Ribose sugar and phosphate
    ]

    if not all(mol.HasSubstructMatch(substructure) for substructure in coa_substructures):
        return False, "Necessary CoA substructures missing"
    
    return True, "Contains trans double bond and key CoA moieties"

# Test examples
smiles_list = [
    "CCCCCCCCCCCCCCCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12",
    "COc1cc(\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2OP(O)(O)=O)n2cnc3c(N)ncnc23)cc(OC)c1O",
]
results = [is_trans_2_enoyl_CoA(smiles) for smiles in smiles_list]
for smiles, (is_match, reason) in zip(smiles_list, results):
    print(f"SMILES: {smiles} is trans-2-enoyl-CoA: {is_match}, Reason: {reason}")