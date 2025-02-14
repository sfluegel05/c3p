"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) based on its SMILES string.
    This class includes acyl-CoA(4-) molecules where the acyl group is a 3-substituted propionyl group,
    and the CoA core is deprotonated at phosphate and diphosphate groups (4- charge).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-substituted propionyl-CoA(4-), False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define key substructures of CoA using SMARTS
    # Adenine moiety
    adenine_smarts = "n1cnc2c(ncnc12)N"
    adenine_pattern = Chem.MolFromSmarts(adenine_smarts)

    # Pantetheine unit with thiol
    pantetheine_smarts = "NCC(=O)NCCS"
    pantetheine_pattern = Chem.MolFromSmarts(pantetheine_smarts)

    # Check for adenine moiety
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Adenine moiety not found (CoA core not present)"

    # Check for pantetheine unit
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Pantetheine unit not found (CoA core not present)"

    # Check for acyl group attached via thioester bond to sulfur
    thioester_smarts = "C(=O)SCCN"
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester bond not found (acyl group not attached to CoA)"

    # Identify acyl chain length and substitution at position 3
    # Find the carbonyl carbon attached to sulfur
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    for match in thioester_matches:
        carbonyl_c_idx = match[0]
        sulfur_idx = match[2]

        # Get the carbonyl carbon atom
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
        # Find the acyl chain connected to the carbonyl carbon (excluding sulfur)
        acyl_chain = []
        stack = [atom for atom in carbonyl_c.GetNeighbors() if atom.GetIdx() != sulfur_idx]
        visited = set()
        while stack:
            atom = stack.pop()
            idx = atom.GetIdx()
            if idx in visited:
                continue
            visited.add(idx)
            acyl_chain.append(atom)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() != carbonyl_c_idx and neighbor.GetIdx() not in visited:
                    stack.append(neighbor)

        # Check if acyl chain is 3 carbons long (propionyl group)
        if len(acyl_chain) < 2:
            continue  # Acyl chain too short
        elif len(acyl_chain) > 2:
            continue  # Acyl chain too long for propionyl group

        # Check for substitution at position 3
        c1 = acyl_chain[0]
        c2 = acyl_chain[1]
        # c2 is the 3rd carbon in the chain (position 3)
        substituents = [nbr for nbr in c2.GetNeighbors() if nbr.GetIdx() != c1.GetIdx()]
        if len(substituents) > 1:
            return True, "Molecule is a 3-substituted propionyl-CoA(4-)"

    return False, "Acyl group is not a 3-substituted propionyl group"