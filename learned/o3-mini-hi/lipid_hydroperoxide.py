"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
"""
Classifies: Lipid Hydroperoxide - Any lipid carrying one or more hydroperoxy substituents.
This heuristic classification checks for the presence of a hydroperoxy (–OOH) group 
and a long aliphatic (lipid-like) chain (at least 8 contiguous carbons).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is defined heuristically as a molecule with at least one hydroperoxy (–OOH) group 
    and a long aliphatic chain (as a proxy for lipid features).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a lipid hydroperoxide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for hydroperoxy group: The hydroperoxy group (–OOH) is represented by the substructure pattern [OX2H]-[OX2]
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2H]-[OX2]")
    if not mol.HasSubstructMatch(hydroperoxy_pattern):
        return False, "No hydroperoxy (–OOH) group found"

    # 2. Check for lipid-like structure: We assume the presence of at least one long aliphatic chain
    # Here we use a simple SMARTS pattern for eight contiguous carbons ("CCCCCCCC")
    lipid_chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    if not mol.HasSubstructMatch(lipid_chain_pattern):
        return False, "No long aliphatic chain found; molecule may not be lipid-like"

    # Optional: Check molecular weight as a rough filter (typical lipids are moderately heavy)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150: # threshold chosen heuristically
        return False, "Molecular weight too low for a typical lipid"

    return True, "Molecule is classified as a lipid hydroperoxide: contains a hydroperoxy group and long aliphatic chain"

# Example usage:
if __name__ == "__main__":
    test_smiles = "O(O)[C@@H](CCCC(O)=O)/C=C/C=C\\C/C=C\\C/C=C\\CC"  # Example: 5S-HpEPE
    result, reason = is_lipid_hydroperoxide(test_smiles)
    print(result, reason)