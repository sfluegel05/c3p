"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: CHEBI:51443 anilide

An anilide is any aromatic amide obtained by acylation of aniline.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the anilide core structure
    anilide_pattern = Chem.MolFromSmarts("c1ccc(cc1)NC(=O)")
    anilide_matches = mol.GetSubstructMatches(anilide_pattern)

    if not anilide_matches:
        return False, "Does not contain the anilide core structure"

    # Ensure the entire aromatic system is aromatic
    aromatic_atoms = [mol.GetAtomWithIdx(idx).GetIsAromatic() for match in anilide_matches for idx in match]
    if not all(aromatic_atoms):
        return False, "The aromatic system is not fully aromatic"

    # Check for common substituents on the aromatic ring(s)
    allowed_substituents = "[Cl,Br,I,F,#17,#35,#53,#16,#7,#6]"  # Halogens, -NH2, -NO2, -OH, alkyls, aryls
    substituent_pattern = Chem.MolFromSmarts(f"c1ccc(cc1)NC(=O)[{allowed_substituents}]")
    substituent_matches = mol.GetSubstructMatches(substituent_pattern)

    if not substituent_matches:
        return False, "The substituents on the aromatic ring(s) are not allowed"

    # Check molecular weight
    mol_weight = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 150 or mol_weight > 500:
        return False, "Molecular weight outside the typical range for anilides"

    return True, "Contains an aromatic amide group connected to an aromatic ring, characteristic of anilides"