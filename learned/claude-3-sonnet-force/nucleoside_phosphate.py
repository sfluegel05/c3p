"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate is a nucleobase-containing molecular entity that is a nucleoside
    in which one or more of the sugar hydroxy groups has been converted into a mono- or poly-phosphate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns for purine and pyrimidine bases
    purine_pattern = Chem.MolFromSmarts("n1cnc2ncnc1[nH]2")
    pyrimidine_pattern = Chem.MolFromSmarts("n1c(=O)[nH]c(N)nc1")

    # Check if molecule contains a nucleobase
    purine_match = mol.HasSubstructMatch(purine_pattern)
    pyrimidine_match = mol.HasSubstructMatch(pyrimidine_pattern)
    if not purine_match and not pyrimidine_match:
        return False, "No nucleobase found"

    # Define pattern for sugar ring
    sugar_pattern = Chem.MolFromSmarts("[OX2r3][CX4r3][CX4r3][CX4r3][CX4r3][OX2r3]")

    # Check if molecule contains a sugar ring
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar ring found"

    # Check if a phosphate group is attached to the sugar ring
    phosphate_attached = False
    for sugar_match in sugar_matches:
        for atom_idx in sugar_match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == "O" and atom.GetIsInRing():
                neighbor_atoms = atom.GetNeighbors()
                for neighbor in neighbor_atoms:
                    if neighbor.GetAtomicNum() == 15:  # Phosphorus
                        phosphate_attached = True
                        break

    if not phosphate_attached:
        return False, "Nucleobase and sugar ring present, but no phosphate group attached"

    # Check for additional phosphate groups (for poly-phosphates)
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) > 1:
        return True, "Contains nucleobase, sugar ring, and poly-phosphate group"
    else:
        return True, "Contains nucleobase, sugar ring, and mono-phosphate group"