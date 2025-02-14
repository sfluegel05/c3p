"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: CHEBI:35621 acyl-CoA
An acyl-CoA is a thioester that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any carboxylic acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for disconnected components
    if len(Chem.GetMolFrags(mol)) > 1:
        components = Chem.GetMolFrags(mol)
        fragments = Chem.FragmentOnBRICSBonds(mol)
        coenzyme_a_fragment = None
        thioester_fragment = None

        for frag in fragments:
            if frag.GetNumAtoms() > 20:
                coenzyme_a_fragment = frag
            else:
                thioester_fragment = frag

        if coenzyme_a_fragment is None or thioester_fragment is None:
            return False, "Cannot identify CoA and thioester components"

        mol = coenzyme_a_fragment

    # Look for coenzyme A moiety
    coenzyme_a_pattern = Chem.MolFromSmarts("[N&R]1C=NC2=C1N=CN=C2N[C@H]3[C@@H]([C@H]([C@@H](O3)COP(=O)(O)O)O)O[C@H]([C@@H]([C@H](OP(=O)(O)O)O)O)COP(=O)(O)O")
    if not mol.HasSubstructMatch(coenzyme_a_pattern):
        return False, "No coenzyme A moiety found"

    # Look for thioester group (-C(=O)-S-)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found"

    # Check if thioester is connected to coenzyme A
    for match in thioester_matches:
        thioester_atom = mol.GetAtomWithIdx(match[1])
        neighbors = [mol.GetAtomWithIdx(nbr).GetSymbol() for nbr in thioester_atom.GetNeighbors()]
        if "N" in neighbors:
            break
    else:
        return False, "Thioester group not connected to coenzyme A moiety"

    # Consider stereochemistry
    Chem.AssignStereochemistryFrom3D(mol)
    Chem.WedgeMolBonds(mol, mol.GetConformer())
    
    # Check molecular weight - acyl-CoAs typically >700 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, "Molecular weight too low for acyl-CoA"

    return True, "Contains coenzyme A moiety with a thioester group attached"