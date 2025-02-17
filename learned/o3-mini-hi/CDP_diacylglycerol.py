"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
"""
Classifies: CDP-diacylglycerol
Definition: A CDP-glycerol having unspecified acyl groups (most commonly fatty acyl groups) at the 1- and 2-positions.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_CDP_diacylglycerol(smiles: str):
    """
    Determines if a molecule is a CDP-diacylglycerol based on its SMILES string.
    A CDP-diacylglycerol is defined as a CDP-glycerol (i.e. a cytidine diphosphate linked to glycerol)
    with acyl groups (usually fatty acyl chains) at the 1- and 2-positions of glycerol.
    
    The function applies heuristic substructure searches:
      1. It looks for a cytidine moiety (assumed to be represented with a nucleobase pattern).
      2. It checks for at least two ester bonds that are candidates for acyl substituents.
      3. It verifies the presence of a diphosphate linkage by finding two phosphorus atoms connected via a bridging oxygen.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a CDP-diacylglycerol, False otherwise.
        str: A message describing the reason for the classification decision.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # 1. Look for a cytidine moiety.
    # This SMARTS looks for a pyrimidine ring with an amino substituent and a carbonyl,
    # which is common in the cytosine part of cytidine.
    cytidine_smarts = "n1ccc(N)nc1=O"
    cytidine_pattern = Chem.MolFromSmarts(cytidine_smarts)
    if not mol.HasSubstructMatch(cytidine_pattern):
        return False, "Cytidine moiety not found"

    # 2. Look for ester bonds (acyl groups). We expect at least two ester linkages on glycerol.
    # The following pattern matches an ester oxygen connected to a carbonyl group.
    ester_smarts = "[OX2H0][CX3](=O)[#6]"  # O-CO-carbon (the acyl part)
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups; at least 2 acyl ester bonds expected"

    # 3. Check for a diphosphate linkage:
    # We need to verify that two phosphorus atoms (atomic number 15) are connected via a bridging oxygen.
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    diphosphate_found = False
    for atom in p_atoms:
        for nbr in atom.GetNeighbors():
            # Look for an oxygen neighbor
            if nbr.GetAtomicNum() == 8:
                # Now check if this oxygen has another P neighbor (other than the original)
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetIdx() != atom.GetIdx() and nbr2.GetAtomicNum() == 15:
                        diphosphate_found = True
                        break
            if diphosphate_found:
                break
        if diphosphate_found:
            break

    if not diphosphate_found:
        return False, "Diphosphate bridging (two P atoms connected via an O) not found"

    # Additional checks (if desired) could inspect the connectivity of the glycerol group.
    # However, because CDP-diacylglycerol structures are complex and heterogeneous,
    # we assume that a cytidine moiety + diphosphate linkage + two esters suffices.

    return True, "Molecule contains a cytidine moiety, diphosphate linkage, and at least 2 acyl ester groups consistent with CDP-diacylglycerol"