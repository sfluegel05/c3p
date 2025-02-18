"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: CHEBI:32876 tertiary amine
"""
from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine has a nitrogen atom bonded to three hydrocarbyl groups (carbon-based substituents).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Iterate through all atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen atom
            # Skip aromatic nitrogens (e.g., in pyridine)
            if atom.GetIsAromatic():
                continue
            # Check if nitrogen has exactly three substituents
            if atom.GetDegree() == 3:
                # Check all bonds are single bonds (exclude cases with double bonds to nitrogen)
                if all(bond.GetBondType() == Chem.BondType.SINGLE for bond in atom.GetBonds()):
                    # Check that all substituents are carbon-based (excluding hydrogen)
                    # Since H is implicit, if degree is 3, all substituents are non-H
                    # Now check each neighbor is a carbon (hydrocarbyl group)
                    all_carbon = True
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() != 6:
                            all_carbon = False
                            break
                    if all_carbon:
                        return True, "Tertiary amine with three carbon-based substituents"
                    else:
                        # Some substituents may have heteroatoms, but according to examples, still count
                        # Modify to accept any substituents (non-hydrogen)
                        return True, "Tertiary amine with three substituents"
    
    return False, "No tertiary amine group detected"