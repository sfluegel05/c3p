"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
"""
Classifies: CHEBI:37412 polyprenol phosphate

A prenol phosphate resulting from the formal condensation of the terminal allylic 
hydroxy group of a polyprenol with 1 mol eq. of phosphoric acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for phosphate group (P(O)(O)=O)
    phosphate_pattern = Chem.MolFromSmarts("P(O)(O)=O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Look for prenol chain (long chain of C=C-C units)
    prenol_pattern = Chem.MolFromSmarts("[CX3]=C[CX3]C=[CX3]")
    prenol_matches = mol.GetSubstructMatches(prenol_pattern)

    # Check for hydroxyl group attached to prenol chain
    has_hydroxy = False
    for match in prenol_matches:
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:  # Oxygen atom
                    has_hydroxy = True
                    break
            if has_hydroxy:
                break

    if not has_hydroxy:
        return False, "No hydroxyl group attached to prenol chain"

    # Count carbon atoms and double bonds to verify chain length
    n_carbon = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)

    if n_carbon < 10 or n_double_bonds < 3:
        return False, "Prenol chain too short"

    # Check molecular weight range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 1000:
        return False, "Molecular weight outside typical range for polyprenol phosphates"

    return True, "Contains a prenol chain with a terminal hydroxyl group condensed with phosphoric acid"