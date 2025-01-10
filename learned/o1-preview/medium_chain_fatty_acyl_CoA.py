"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: medium-chain fatty acyl-CoA
Definition: A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any medium-chain fatty acid.
"""

from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify the thioester bond: C(=O)-S
    pattern = Chem.MolFromSmarts("C(=O)S")
    matches = mol.GetSubstructMatches(pattern)
    if len(matches) == 0:
        return False, "Thioester linkage to CoA not found"

    # Iterate over all thioester matches
    for match in matches:
        c_idx = match[0]  # Carbonyl carbon atom index
        s_idx = match[2]  # Sulfur atom index

        # Get the bond between carbonyl carbon and sulfur
        bond = mol.GetBondBetweenAtoms(c_idx, s_idx)
        if bond is None:
            continue  # Skip if bond is not found

        # Break the bond between carbonyl carbon and sulfur
        bond_idx = bond.GetIdx()
        fragments = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=True)
        frags = Chem.GetMolFrags(fragments, asMols=True, sanitizeFrags=True)

        # Identify the acyl chain fragment and CoA fragment
        acyl_chain = None
        coa_moiety = None

        for frag in frags:
            atom_nums = {atom.GetAtomicNum() for atom in frag.GetAtoms()}
            # CoA contains phosphorus (15) and multiple nitrogen atoms (7)
            if 15 in atom_nums and atom_nums.count(7) >= 5:
                coa_moiety = frag
            elif 6 in atom_nums and 15 not in atom_nums:
                acyl_chain = frag

        if acyl_chain is None or coa_moiety is None:
            continue  # Try next thioester linkage

        # Count the number of carbons in the acyl chain
        num_carbons = sum(1 for atom in acyl_chain.GetAtoms() if atom.GetAtomicNum() == 6)

        if num_carbons < 6:
            return False, f"Acyl chain too short ({num_carbons} carbons), not medium-chain"
        if num_carbons > 12:
            return False, f"Acyl chain too long ({num_carbons} carbons), not medium-chain"

        # Check for key features of CoA in coa_moiety

        # Adenine ring pattern
        adenine_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)")
        if not coa_moiety.HasSubstructMatch(adenine_pattern):
            return False, "Adenine moiety not found in CoA fragment"

        # Diphosphate linkage pattern
        diphosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)OP(=O)(O)O")
        if not coa_moiety.HasSubstructMatch(diphosphate_pattern):
            return False, "Diphosphate linkage not found in CoA fragment"

        # Pantetheine moiety pattern
        pantetheine_pattern = Chem.MolFromSmarts("SCCNC(=O)CC(O)C(C)(C)")
        if not coa_moiety.HasSubstructMatch(pantetheine_pattern):
            return False, "Pantetheine moiety not found in CoA fragment"

        # If all checks pass, classify as medium-chain fatty acyl-CoA
        return True, f"Contains medium-chain acyl group ({num_carbons} carbons) attached to Coenzyme A via thioester linkage"

    # If no valid thioester linkage is found
    return False, "Does not contain medium-chain fatty acyl-CoA moiety"