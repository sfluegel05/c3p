"""
Classifies: CHEBI:140310 phenyl acetates
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is an acetate ester obtained by formal condensation
    of the carboxy group of acetic acid with the hydroxy group of any phenol.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a phenyl acetate, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define acetate ester part
    # Represents ester linkage with aryl group
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    
    # Define a generic phenyl (aryl) pattern with potential for attachment
    phenol_derived_pattern = Chem.MolFromSmarts("c[OX2H0]")

    # Check if the molecule has an ester pattern
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"
    
    # Check if the molecule has an aryl structure with a hydroxy connection
    if not mol.HasSubstructMatch(phenol_derived_pattern):
        return False, "No aryl structure with hydroxy connection found"

    # Ensure both structures are present within the same molecule
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    aryl_matches = mol.GetSubstructMatches(phenol_derived_pattern)

    for ester in ester_matches:
        ester_ester_oxygen = ester[0]  # oxygen in ester

        for aryl in aryl_matches:
            aryl_hydroxy_oxygen = aryl[1]  # oxy connection in aryl structure

            # Check if there is a bond between the ester oxygen and aryl hydroxy oxygen
            if mol.GetBondBetweenAtoms(ester_ester_oxygen, aryl_hydroxy_oxygen) is not None:
                return True, "Contains a linkage of an acetate ester to an aryloxy group"

    return False, "Ester linkage and aryl structure are not properly connected"