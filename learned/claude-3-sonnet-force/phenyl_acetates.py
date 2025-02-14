"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: CHEBI:37405 phenyl acetates
An acetate ester obtained by formal condensation of the carboxy group of acetic acid
with the hydroxy group of any phenol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenyl acetate, False otherwise
        str: Reason for classification
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Look for aromatic ring with hydroxy group
        phenol_pattern = Chem.MolFromSmarts("c1ccc(cc1)O")
        phenol_matches = mol.GetSubstructMatches(phenol_pattern)
        if not phenol_matches:
            return False, "No phenol substructure found"

        # Look for acetate group (-O-C(=O)-C)
        acetate_pattern = Chem.MolFromSmarts("OC(=O)C")
        acetate_matches = mol.GetSubstructMatches(acetate_pattern)
        if not acetate_matches:
            return False, "No acetate group found"

        # Check if acetate group is attached to phenol ring
        for acetate_match in acetate_matches:
            acetate_oxygen = mol.GetAtomWithIdx(acetate_match[0])
            for neighbor in acetate_oxygen.GetNeighbors():
                if neighbor.IsInRing() and neighbor.IsAromatic():
                    # Additional checks for phenyl acetate definition
                    acetic_acid_pattern = Chem.MolFromSmarts("CC(=O)O")
                    if mol.HasSubstructMatch(acetic_acid_pattern):
                        return True, "Contains a phenol ring with an acetate group attached"
                    else:
                        return False, "Acetate group not derived from acetic acid"

        return False, "Acetate group not attached to phenol ring"

    except Exception as e:
        return False, f"Error occurred: {str(e)}"