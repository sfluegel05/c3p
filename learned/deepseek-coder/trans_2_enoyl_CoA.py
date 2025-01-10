"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: CHEBI:28348 trans-2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    A trans-2-enoyl-CoA is characterized by a CoA moiety and a trans-2-enoyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trans-2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible CoA moiety pattern
    coa_pattern = Chem.MolFromSmarts(
        "[*]SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
    )
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Define the trans-2-enoyl pattern with more specificity
    # - Must have exactly one trans double bond at position 2
    # - Must have exactly one carbonyl at position 1
    # - Must be connected to CoA via thioester bond
    trans_2_enoyl_pattern = Chem.MolFromSmarts("[CX3H1,CX4H2]\\C=C\\C(=O)[SX2]")
    matches = mol.GetSubstructMatches(trans_2_enoyl_pattern)
    
    if not matches:
        return False, "No trans-2-enoyl group found"
    
    # Verify there's exactly one trans-2-enoyl group
    if len(matches) != 1:
        return False, f"Found {len(matches)} potential enoyl groups, need exactly 1"

    # Verify the double bond is truly trans
    bond = mol.GetBondBetweenAtoms(matches[0][1], matches[0][2])
    if not bond:
        return False, "No double bond found between expected atoms"
    
    # Check if stereo information is present and correct
    if bond.GetStereo() == Chem.rdchem.BondStereo.STEREOZ:
        return False, "Double bond is in cis configuration"
    
    # If no stereo information, check the geometry
    if bond.GetStereo() == Chem.rdchem.BondStereo.STEREONONE:
        # Use 3D coordinates to determine trans configuration
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        conf = mol.GetConformer()
        coords1 = conf.GetAtomPosition(matches[0][1])
        coords2 = conf.GetAtomPosition(matches[0][2])
        coords3 = conf.GetAtomPosition(matches[0][3])
        angle = Chem.rdMolTransforms.GetAngleRad(conf, matches[0][1], matches[0][2], matches[0][3])
        if abs(angle) < 2.8:  # Approximate threshold for trans configuration
            return False, "Double bond geometry suggests cis configuration"

    # Verify the thioester bond is connected to CoA
    thioester_atom = matches[0][-1]  # The sulfur atom
    for neighbor in mol.GetAtomWithIdx(thioester_atom).GetNeighbors():
        if neighbor.GetSymbol() == "C" and neighbor.GetDegree() == 3:  # CoA connection
            return True, "Contains CoA moiety and trans-2-enoyl group with a thioester bond"

    return False, "Thioester bond not properly connected to CoA moiety"