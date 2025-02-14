"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:60903 medium-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.
    A medium-chain fatty acyl-CoA results from the condensation of coenzyme A with a medium-chain fatty acid.
    Medium-chain fatty acids are aliphatic chains with 6 to 12 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the Coenzyme A (CoA) SMARTS pattern
    # This pattern captures the core structure of CoA, including the phosphates, ribose, adenine, and pantetheine moiety
    coa_smarts = Chem.MolFromSmarts('NC(=O)CCNC(=O)[C@@H](O)C(C)(C)[C@@H](OP(=O)(O)OP(=O)(O)OC[C@@H]1O[C@H](CO[P](=O)(O)O)[C@H](O)[C@@H]1O)O')

    # Check if the molecule contains the CoA moiety
    if not mol.HasSubstructMatch(coa_smarts):
        return False, "Coenzyme A moiety not found"

    # Define the thioester linkage pattern connecting the acyl chain to CoA
    # We look for the C(=O)-S linkage where S is connected to CoA
    thioester_pattern = Chem.MolFromSmarts('C(=O)SCCNC(=O)')

    # Find the thioester linkage
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester linkage to CoA not found"

    # Assume the first match is the correct thioester linkage
    thioester_match = thioester_matches[0]
    carbonyl_c_idx = thioester_match[0]
    sulfur_idx = thioester_match[1]

    # Break the bond between the carbonyl carbon and the sulfur atom to isolate the acyl chain
    emol = Chem.EditableMol(mol)
    emol.RemoveBond(carbonyl_c_idx, sulfur_idx)
    mol_frag = emol.GetMol()

    # Get the fragments of the molecule after the bond break
    fragments = Chem.GetMolFrags(mol_frag, asMols=True)

    # Identify the acyl chain fragment by searching for the fragment containing the carbonyl carbon
    acyl_chain = None
    for frag in fragments:
        # The acyl chain should contain the carbonyl group C=O
        if frag.HasSubstructMatch(Chem.MolFromSmarts('C=O')):
            acyl_chain = frag
            break

    if acyl_chain is None:
        return False, "Acyl chain not found"

    # Count the number of carbon atoms in the acyl chain, including the carbonyl carbon
    carbon_count = 0
    for atom in acyl_chain.GetAtoms():
        if atom.GetAtomicNum() == 6:
            carbon_count += 1

    # Medium-chain fatty acids have 6 to 12 carbons
    if 6 <= carbon_count <= 12:
        return True, f"Contains fatty acyl chain with {carbon_count} carbons (medium-chain)"
    else:
        return False, f"Fatty acyl chain with {carbon_count} carbons is not medium-chain"