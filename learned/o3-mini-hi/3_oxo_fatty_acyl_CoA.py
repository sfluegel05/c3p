"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: 3-oxo-fatty acyl-CoA
An oxo fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A 
with the carboxy group of any 3-oxo-fatty acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    The defined criteria are:
      1. Presence of a 3-oxo fatty acyl thioester moiety. This is captured by the substructure
         pattern [CX3](=O)[CH2]C(=O)S representing R-C(=O)-CH2-C(=O)-S.
      2. Presence of a coenzyme A moiety. One of the hallmarks is the adenine substructure 
         (as part of the ADP portion of CoA). We capture a typical adenine using the SMARTS pattern
         n1cnc2ncnc(c12).
    
    Args:
        smiles (str): A SMILES string representing the molecule.
    
    Returns:
        bool: True if the molecule is a 3-oxo-fatty acyl-CoA, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the substructure pattern for the 3-oxo fatty acyl thioester fragment.
    # This pattern looks for a carbonyl group directly attached to a CH2 that is connected to 
    # a second carbonyl which is bonded to a sulfur atom.
    acyl_pattern = Chem.MolFromSmarts("[CX3](=O)[CH2]C(=O)S")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "Does not contain a 3-oxo fatty acyl thioester substructure"

    # Define a substructure pattern for a key adenine fragment that is part of Coenzyme A.
    # Many 3-oxo-fatty acyl-CoA molecules contain this adenine ring.
    coa_pattern = Chem.MolFromSmarts("n1cnc2ncnc(c12)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Does not contain a Coenzyme A adenine fragment"

    # Optionally one can further check for phosphate groups or other CoA-specific fragments.
    # For the purposes of classification, the above two matches are sufficient.
    return True, "Contains a 3-oxo fatty acyl thioester linked to a coenzyme A moiety"