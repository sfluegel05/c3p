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
    
    # Parse SMILES to RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Attempt to define SMARTS patterns for trans double bond and CoA components
    try:
        trans_double_bond_patterns = [
            Chem.MolFromSmarts("C/C=C\\C"),  # Trans-double bonds
            Chem.MolFromSmarts("C\\C=C/C")
        ]
        
        coa_substructures = [
            Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)C"),  # Key CoA linkage and structure
            Chem.MolFromSmarts("COP(O)(=O)OP(O)(=O)O"),  # Phosphate groups
            Chem.MolFromSmarts("OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")  # Ribose and nucleotide
        ]
        
        # Validate substructure patterns
        for pattern in trans_double_bond_patterns + coa_substructures:
            if pattern is None:
                raise ValueError("Invalid SMARTS pattern initialization")
    
        # Match trans double bond
        if not any(mol.HasSubstructMatch(pattern) for pattern in trans_double_bond_patterns):
            return False, "No trans double bond found"
        
        # Match CoA components
        if not all(mol.HasSubstructMatch(substructure) for substructure in coa_substructures):
            return False, "Necessary CoA substructures missing"
    
        return True, "Contains trans double bond and key CoA moieties"

    except Exception as e:
        return False, f"Error during matching: {str(e)}"

# Test examples with logging of results
smiles_list = [
    "CCCCCCCCCCCCCCCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(=O)=O)n1cnc2c(N)ncnc12",
    "COc1cc(\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2OP(O)(O)=O)n2cnc3c(N)ncnc23)cc(OC)c1O",
]

results = [is_trans_2_enoyl_CoA(smiles) for smiles in smiles_list]
for smiles, (is_match, reason) in zip(smiles_list, results):
    print(f"SMILES: {smiles} is trans-2-enoyl-CoA: {is_match}, Reason: {reason}")