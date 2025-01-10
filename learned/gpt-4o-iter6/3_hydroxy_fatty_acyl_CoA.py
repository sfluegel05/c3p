"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA is characterized by a 3-hydroxy fatty acid linked to coenzyme A.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Update SMARTS pattern for Coenzyme A - focus on key motifs (phosphate, adenosine, pantetheine)
    coA_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)C1CNC2=C1N=CN=C2N3[C@H]4O[C@@H](COP(O)(=O)OP(O)(O)=O)[C@H]4O[C@H]3CO")
    
    # Updated generalized SMARTS pattern for 3-hydroxy fatty acid: flexible chain with C=O and C-OH separation
    hydroxy_fatty_acid_pattern = Chem.MolFromSmarts("C[C@@H](O)[CX3](=O)*)")

    # Check for Coenzyme A moiety
    if not mol.HasSubstructMatch(coA_pattern):
        return False, "No Coenzyme A moiety found"

    # Check for 3-hydroxy fatty acid component
    if not mol.HasSubstructMatch(hydroxy_fatty_acid_pattern):
        return False, "No 3-hydroxy fatty acid moiety found"

    return True, "Molecule contains both Coenzyme A and 3-hydroxy fatty acid moieties"