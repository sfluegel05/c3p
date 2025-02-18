"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: CHEBI:3616 furanocoumarin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin consists of a furan ring fused with a coumarin (benzene fused to a pyrone ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define coumarin substructure (benzene fused to pyrone: c1ccc2oc(=O)cc2c1)
    coumarin_pattern = Chem.MolFromSmarts('c1ccc2oc(=O)cc2c1')
    coumarin_matches = mol.GetSubstructMatches(coumarin_pattern)
    if not coumarin_matches:
        return False, "No coumarin substructure found"

    # Define furan substructure (o1cccc1 with aromatic bonds)
    furan_pattern = Chem.MolFromSmarts('o1:c:c:c:c1')
    furan_matches = mol.GetSubstructMatches(furan_pattern)
    if not furan_matches:
        return False, "No furan ring found"

    # Check fusion between any coumarin and furan rings
    for cm in coumarin_matches:
        coumarin_atoms = set(cm)
        for fm in furan_matches:
            furan_atoms = set(fm)
            shared = coumarin_atoms & furan_atoms
            # Fusion requires at least two shared atoms (adjacent in both rings)
            if len(shared) >= 2:
                # Verify adjacency in the furan ring
                furan_bonds = mol.GetBondsBetweenAtoms(fm[0], fm[1])
                for i in range(len(fm)-1):
                    a1 = fm[i]
                    a2 = fm[i+1]
                    if a1 in shared and a2 in shared:
                        return True, "Furan ring fused to coumarin structure"
                # Check last to first atom in the ring
                a1 = fm[-1]
                a2 = fm[0]
                if a1 in shared and a2 in shared:
                    return True, "Furan ring fused to coumarin structure"

    return False, "Furan not fused to coumarin"