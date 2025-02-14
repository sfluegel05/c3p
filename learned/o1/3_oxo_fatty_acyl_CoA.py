"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
"""
Classifies: 3-oxo-fatty acyl-CoA
Definition: An oxo fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any 3-oxo-fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-oxo-fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Remove hydrogen atoms for substructure matching
    mol = Chem.RemoveHs(mol)

    # Define the CoA substructure pattern (more comprehensive)
    # Using the SMILES of Coenzyme A minus the thiol hydrogen
    coa_smiles = "CC(C)(COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n2cnc3c(N)ncnc23)O[C@@H](C(=O)NCCC(=O)NCCS)C(=O)NCCS"
    coa_mol = Chem.MolFromSmiles(coa_smiles)
    if not coa_mol:
        return False, "Invalid CoA SMILES"

    # Check for CoA substructure
    if not mol.HasSubstructMatch(coa_mol):
        return False, "CoA moiety not found"

    # Define the thioester linkage pattern
    thioester_smarts = "C(=O)SCCNC(=O)"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage not found"

    # Define the 3-oxo-fatty acyl chain pattern
    # Pattern for thioester-linked 3-oxo-fatty acyl chain
    # The pattern starts from the carbonyl carbon of the thioester linkage
    # and includes the beta-keto group at the 3-position
    oxo_fatty_acyl_smarts = "C(=O)SC[C;H2][C](=O)[C;H]"
    oxo_fatty_acyl_pattern = Chem.MolFromSmarts(oxo_fatty_acyl_smarts)
    if not mol.HasSubstructMatch(oxo_fatty_acyl_pattern):
        return False, "3-oxo-fatty acyl chain not found"

    # Optional: Check the length of the fatty acyl chain
    # Exclude very short chains that are not considered fatty acids
    chain_length = 0
    for match in mol.GetSubstructMatches(oxo_fatty_acyl_pattern):
        # The match includes the thioester linkage and the beta-keto group
        # Count the number of carbons in the chain beyond the beta carbon
        beta_carbon_idx = match[3]
        chain_atoms = Chem.rdmolops.GetShortestPath(mol, beta_carbon_idx, beta_carbon_idx)
        for atom_idx in chain_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:
                chain_length += 1
        break  # Only need the first match

    if chain_length < 4:
        return False, f"Fatty acyl chain too short ({chain_length} carbons)"

    return True, "Molecule is a 3-oxo-fatty acyl-CoA"