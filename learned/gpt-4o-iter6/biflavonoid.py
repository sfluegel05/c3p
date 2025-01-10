"""
Classifies: CHEBI:50128 biflavonoid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    A biflavonoid is defined as a flavonoid oligomer obtained by the oxidative coupling 
    of at least two units of aryl-substituted benzopyran rings or its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a biflavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for flavonoid and connection
    flavonoid_unit_pattern = Chem.MolFromSmarts("c1c(O)c2cc(O)ccc2oc1")  # captures the benzopyran framework
    connection_pattern = Chem.MolFromSmarts("c-c")
    
    # Match flavonoid units
    flavonoid_matches = mol.GetSubstructMatches(flavonoid_unit_pattern)
    if len(flavonoid_matches) < 2:
        return False, "Less than two flavonoid units found"
    
    # Check for direct connectivity between flavonoid units
    for i in range(len(flavonoid_matches)):
        unit1_atoms = {atom_idx for atom_idx in flavonoid_matches[i]}
        for j in range(i+1, len(flavonoid_matches)):
            unit2_atoms = {atom_idx for atom_idx in flavonoid_matches[j]}
            # Check if there's a bond connecting these two sets
            connected = any(
                mol.GetBondBetweenAtoms(idx1, idx2)
                for idx1 in unit1_atoms
                for idx2 in unit2_atoms
                if mol.GetBondBetweenAtoms(idx1, idx2) is not None
            )
            if connected:
                return True, "Contains at least two flavonoid units linked by a single bond/atom"

    return False, "No connecting bond/atom found between flavonoid units"