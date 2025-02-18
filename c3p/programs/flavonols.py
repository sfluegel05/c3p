"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: flavonols
Definition: Any hydroxyflavone in which the ring hydrogen at position 3 of the heterocyclic ring is replaced by a hydroxy group.
In other words, the molecule must contain a 3-hydroxyflavone (flavonol) core.
Examples include tambulin, pinoquercetin, etc.
"""

from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    Flavonols are defined as hydroxyflavones in which the ring hydrogen at position 3
    of the heterocyclic (chromen-4-one) ring is replaced by a hydroxy group.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a flavonol, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved SMARTS pattern for the 3-hydroxyflavone core.
    # This pattern looks for:
    #   - A benzene ring (ring A) fused to a heterocycle (ring C) via the pattern "c1ccc2c(c1)"
    #   - The heterocycle must contain:
    #         • an oxygen atom as part of the ring ("oc")
    #         • a carbonyl group ("c2=O") and 
    #         • a carbon that carries a hydroxy group ([OX2H]) in place of the ring hydrogen, i.e. at position 3.
    # The resulting SMARTS is:
    flavonol_pattern = Chem.MolFromSmarts("c1ccc2c(c1)oc(c([OX2H])c2=O)")
    if flavonol_pattern is None:
        return False, "Error in generating SMARTS pattern"

    # Check if the molecule contains the flavonol scaffold
    if mol.HasSubstructMatch(flavonol_pattern):
        return True, "Molecule contains a 3-hydroxyflavone core (flavonol scaffold)"
    else:
        return False, "Molecule does not contain a 3-hydroxyflavone core required for flavonols"

# Example usage (for testing when run directly):
if __name__ == "__main__":
    # Test example: tambulin (should be classified as flavonol)
    test_smiles = "COc1ccc(cc1)-c1oc2c(OC)c(OC)cc(O)c2c(=O)c1O"
    classification, explanation = is_flavonols(test_smiles)
    print(f"SMILES: {test_smiles}")
    print("Classification:", classification)
    print("Reason:", explanation)