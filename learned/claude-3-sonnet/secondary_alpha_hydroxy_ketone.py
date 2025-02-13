"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: CHEBI:33762 secondary alpha-hydroxy ketone
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone has a carbonyl group and a hydroxy group linked
    by a carbon bearing one hydrogen and one organyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Find all carbonyl and hydroxy groups
        carbonyl_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0]
        hydroxy_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0 and atom.GetTotalNumHs() == 1]

        # Check if any carbonyl and hydroxy groups are adjacent
        for c_idx in carbonyl_atoms:
            for h_idx in hydroxy_atoms:
                if mol.GetBondBetweenAtoms(c_idx, h_idx) is not None:
                    # Check if the shared carbon has one hydrogen and one organyl group
                    shared_carbon_idx = mol.GetBondBetweenAtoms(c_idx, h_idx).GetBeginAtomIdx()
                    if mol.GetAtomWithIdx(shared_carbon_idx).GetHybridization() == Chem.HybridizationType.SP3 and \
                       mol.GetAtomWithIdx(shared_carbon_idx).GetTotalNumHs() == 1 and \
                       len([bond.GetIdx() for bond in mol.GetAtomWithIdx(shared_carbon_idx).GetBonds() if bond.GetBondTypeAsDouble() > 1]) == 1:
                        return True, "Contains a secondary alpha-hydroxy ketone moiety"

        return False, "No secondary alpha-hydroxy ketone moiety found"

    except Exception as e:
        return None, f"An error occurred: {str(e)}"