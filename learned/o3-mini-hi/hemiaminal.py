"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: Hemi-aminal compounds.
Definition: A hemiaminal is any organic compound containing a tetrahedral (sp3, non‐aromatic) carbon with exactly one hydrogen,
that has exactly three heavy-atom neighbors, two of which are substituents: one hydroxyl group (-OH) and one amino group (-NH, -NHR, or -NR2).
This pattern is captured using a SMARTS query.
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule contains a hemiaminal motif.

    A hemiaminal is defined as a tetrahedral (sp3, non-aromatic) carbon carrying exactly one hydrogen and three heavy-atom neighbors,
    with two required substituents: one hydroxyl group (-OH) and one amino group (-NH, -NHR, or -NR2).
    
    This function uses a SMARTS substructure query to locate such a motif.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the hemiaminal motif is detected, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are correct.
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for a hemiaminal motif:
    #   [C;H1;X4]     : a tetrahedral (sp3) carbon with exactly one hydrogen.
    #   ([OX2H])      : one substituent is an oxygen that is bound as hydroxyl (-OH).
    #   ([NX3])       : one substituent is a trivalent nitrogen (covers primary, secondary or tertiary amino groups).
    #   (*)            : the remaining substituent (any heavy atom).
    #
    # This SMARTS enforces that the candidate carbon has exactly three heavy-atom neighbors.
    hemiaminal_smarts = "[C;H1;X4]([OX2H])([NX3])(*)"
    pattern = Chem.MolFromSmarts(hemiaminal_smarts)
    if pattern is None:
        return False, "Internal error: Invalid SMARTS pattern"
    
    # Search for the hemiaminal motif in the molecule.
    if mol.HasSubstructMatch(pattern):
        return True, ("Molecule contains a hemiaminal motif: a tetrahedral sp3 carbon carrying exactly one hydrogen "
                      "with an –OH substituent and an amino substituent.")
    else:
        return False, ("No hemiaminal motif (tetrahedral sp3 carbon with one hydrogen, -OH and -NH substituents) was found.")

# Example usage (uncomment to test):
# test_smiles_list = [
#     "NC(O)C(O)=O",         # alpha-hydroxyglycine (should be True)
#     "OC(N)CC",             # 2-Aminopropanol (should be True)
#     "C1CCCCC1"             # Cyclohexane (should be False)
# ]
# for smi in test_smiles_list:
#     result, explanation = is_hemiaminal(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nExplanation: {explanation}\n")