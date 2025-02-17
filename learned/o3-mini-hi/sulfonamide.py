"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: Sulfonamide 
Definition: An amide of a sulfonic acid RS(=O)2NR'2.
Improved pattern demands:
  – A sulfonyl sulfur (S) that is tetravalent (X4) and double‐bonded to two oxygens
  – At least one carbon substituent on S (i.e. R–S(=O)(=O))
  – An acyclic nitrogen attached to S (i.e. –N not part of a ring), representing the “amide” N.
This extra structure should help avoid matches in molecules where the S(=O)(=O)–N fragment 
appears in an inappropriate context.
"""

from rdkit import Chem

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide is defined as: RS(=O)2NR'2, meaning the sulfur has two double bonds to oxygen,
    is linked to an R group (typically a carbon) and an N that is not part of a ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule contains a sulfonamide group, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # An improved SMARTS that looks for:
    # 1. A carbon group attached to a tetravalent sulfur having two double-bonded oxygens.
    # 2. The sulfur is connected to a nitrogen that is not in a ring.
    # This enforces the typical RS(=O)(=O)-N connectivity.
    improved_pattern = Chem.MolFromSmarts("[#6]-[SX4](=[O])(=[O])-[NX3;!R]")
    if improved_pattern is None:
        return False, "Error in generating SMARTS pattern"
    
    # Find all substructure matches.
    matches = mol.GetSubstructMatches(improved_pattern)
    if matches:
        return True, "Molecule contains a sulfonamide group: RS(=O)(=O)-N (with N acyclic and S attached to a carbon)"
    else:
        return False, "Molecule does not contain the required sulfonamide group RS(=O)(=O)-N (with proper connectivity)"

# If needed, one might add a main block for simple testing:
if __name__ == "__main__":
    # Test example sulfonamide SMILES (from the true positives above)
    test_smiles = "CN(CCOc1ccc(NS(C)(=O)=O)cc1)CCc1ccc(NS(C)(=O)=O)cc1"  # dofetilide
    result, reason = is_sulfonamide(test_smiles)
    print("Result:", result)
    print("Reason:", reason)