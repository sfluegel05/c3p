"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: acyl-CoA
"""
from rdkit import Chem

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is a thioester resulting from the condensation of the thiol group of coenzyme A
    with the carboxy group of any carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES of input molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Coenzyme A SMILES string
    coa_smiles = "CC(=O)NCCSCCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3nc(N)nc(N)c32)C(O)C1OP(=O)(O)O"
    coa_mol = Chem.MolFromSmiles(coa_smiles)
    if coa_mol is None:
        return False, "Error in loading CoA structure"

    # Find sulfur atom index in CoA molecule
    sulfur_idx = None
    for atom in coa_mol.GetAtoms():
        if atom.GetSymbol() == 'S':
            sulfur_idx = atom.GetIdx()
            break
    if sulfur_idx is None:
        return False, "No sulfur atom found in CoA structure"

    # Search for CoA substructure in input molecule
    matches = mol.GetSubstructMatches(coa_mol, useChirality=True)
    if not matches:
        return False, "Coenzyme A moiety not found"

    # For each match, check for thioester linkage
    for match in matches:
        # Get the corresponding sulfur atom in the input molecule
        mol_sulfur_idx = match[sulfur_idx]
        sulfur_atom = mol.GetAtomWithIdx(mol_sulfur_idx)

        # Check neighbors of sulfur atom
        neighbors = sulfur_atom.GetNeighbors()
        found_thioester = False
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 6:  # Carbon
                # Check if this carbon is a carbonyl carbon (C=O)
                is_carbonyl = False
                for nbr_of_carb in neighbor.GetNeighbors():
                    if nbr_of_carb.GetAtomicNum() == 8:  # Oxygen
                        bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nbr_of_carb.GetIdx())
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            is_carbonyl = True
                            break
                if is_carbonyl:
                    found_thioester = True
                    break
        if found_thioester:
            return True, "Contains CoA moiety attached via thioester linkage to an acyl group"

    return False, "Thioester linkage not connected to CoA sulfur atom"