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
      1. Presence of a 3-oxo fatty acyl thioester moiety. Previously the pattern required an exact
         "[CX3](=O)[CH2]C(=O)S" (as in acetoacetyl-CoA), but that misses substituted or cyclic systems.
         Here we use a more flexible pattern "[CX3](=O)[#6]C(=O)S" so that the CH2 can be replaced by any aliphatic carbon.
      2. Presence of a Coenzyme A moiety. A key hallmark is the adenine fragment. We capture a typical adenine
         using the SMARTS pattern "n1cnc2ncnc(c12)".
      3. To avoid classifying ionized (deprotonated) CoA derivatives, we check that the SMILES does not contain any
         explicit deprotonated oxygens "[O-]".
         
    Args:
        smiles (str): SMILES string representing the molecule.
    
    Returns:
        bool: True if the molecule matches the 3-oxo fatty acyl-CoA criteria, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude molecules that have explicit deprotonated phosphate groups (which appear as [O-]) 
    # since we expect a neutral CoA moiety.
    if "[O-]" in smiles:
        return False, "Contains deprotonated phosphate groups, not a neutral CoA moiety"
    
    # Define a substructure SMARTS for the 3-oxo fatty acyl thioester fragment.
    # This pattern now allows any aliphatic carbon (using [#6]) between the two carbonyls.
    acyl_pattern = Chem.MolFromSmarts("[CX3](=O)[#6]C(=O)S")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "Does not contain a 3-oxo fatty acyl thioester substructure"
    
    # Define the adenine substructure pattern that is part of the coenzyme A moiety.
    coa_pattern = Chem.MolFromSmarts("n1cnc2ncnc(c12)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Does not contain a Coenzyme A adenine fragment"
    
    return True, "Contains a 3-oxo fatty acyl thioester linked to a coenzyme A moiety"

# Below is an example of usage:
if __name__ == "__main__":
    # Example SMILES for 3-oxoadipyl-CoA:
    smiles_example = "CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)CC(=O)CCC(O)=O"
    result, reason = is_3_oxo_fatty_acyl_CoA(smiles_example)
    print(result, reason)