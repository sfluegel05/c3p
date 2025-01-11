"""
Classifies: CHEBI:36498 galactosylceramide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide is a cerebroside that contains a galactose head group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a galactosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Recognize galactose via typical 5-ring structure with OH groups
    galactose_pattern = Chem.MolFromSmarts("OC1C(O)C(O)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No galactose residue found"

    # Recognize amide linkage pattern (N-C(=O))
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"

    # Count long hydrocarbon chains which are hallmarks of ceramides
    n_carbon_chains = rdMolDescriptors.CalcNumAliphaticCarbocycles(mol)
    if n_carbon_chains == 0:
        return False, "Long hydrocarbon chains of ceramides not found"

    # Check for presence of a sphingosine or similar long-chain base component
    sphingosine_pattern = Chem.MolFromSmarts("[C;R0][C;R0][C;R0]", mergeHs=True)
    sphingosine_matches = mol.GetSubstructMatches(sphingosine_pattern)
    if len(sphingosine_matches) < 1:
        return False, "Sphingosine-like backbone not detected"

    return True, "Contains galactose, amide linkage, and sphingosine backbone"


__metadata__ = {   'chemical_class': {   'name': 'galactosylceramide',
                          'definition': 'Any of the cerebrosides in which the '
                                        'monosaccharide head group is galactose.',
                          'examples': [
                              "CCCCCCCCCCCCCCCCC(=O)N[C@@H](COC1O[C@H](CO) [C@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC",
                              "CCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COC1O[C@H](CO) [C@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
                          ]
    }
}