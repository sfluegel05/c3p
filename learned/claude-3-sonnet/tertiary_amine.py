"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: CHEBI:33899 tertiary amine
A compound formally derived from ammonia by replacing three hydrogen atoms by hydrocarbyl groups.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondType

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get nitrogen atoms
    n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]

    for n_atom in n_atoms:
        # Exclude quaternary nitrogen and N-oxides
        if n_atom.GetTotalDegree() != 3 or n_atom.GetFormalCharge() != 0:
            continue

        # Count carbon substituents
        c_substituents = 0
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                aromatic = neighbor.GetIsAromatic()
                if aromatic:
                    c_substituents += 1
                else:
                    # Check if the substituent is a hydrocarbyl group
                    for atom in mol.GetAtomWithIdx(neighbor.GetIdx()).GetNeighbors():
                        if atom.GetAtomicNum() != 1 and atom.GetAtomicNum() != 6:
                            break
                    else:
                        c_substituents += 1

        # Check if it has exactly 3 carbon substituents
        if c_substituents == 3:
            return True, "Contains a nitrogen atom with three hydrocarbyl substituents"

    return False, "No tertiary amine group found"