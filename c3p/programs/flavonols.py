"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: flavonols
Definition: Any hydroxyflavone in which the ring hydrogen at position 3 of the heterocyclic ring is replaced by a hydroxy group.
In other words, the molecule must contain a 3-hydroxyflavone core.
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
        bool: True if molecule is a flavonol, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the 3-hydroxyflavone core.
    # Explanation:
    # "c1ccc2c(c1)"   : a benzene ring (ring A) fused to ring C
    # "oc(c(O)c2=O)"   : the heterocyclic ring (ring C) containing:
    #                    - an oxygen in the ring,
    #                    - a carbon bearing a hydroxy group (the 3-OH) and
    #                    - a carbonyl group (at position 4, i.e. "=O")
    #
    # This pattern is designed to match the scaffold when a hydrogen at position 3 is replaced by an –OH.
    flavonol_pattern = Chem.MolFromSmarts("c1ccc2c(c1)oc(c(O)c2=O)")
    if flavonol_pattern is None:
        return False, "Error in generating SMARTS pattern"
    
    if mol.HasSubstructMatch(flavonol_pattern):
        return True, "Molecule contains a 3-hydroxyflavone core (flavonol scaffold)"
    else:
        return False, "Molecule does not contain a 3-hydroxyflavone core required for flavonols"

# Example usage:
if __name__ == "__main__":
    # Test with tambulin SMILES example from the prompt
    test_smiles = "COc1ccc(cc1)-c1oc2c(OC)c(OC)cc(O)c2c(=O)c1O"
    is_class, reason = is_flavonols(test_smiles)
    print(is_class, reason)